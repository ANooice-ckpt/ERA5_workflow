ERA5 server bundle — memory-aware scheduler v2.1
=============================================
更新日期：2026-10-14。完整 Windows→Ubuntu 部署、正式运行和恢复命令见 DEPLOY_MEMORY.txt。
此版本复用现有 models/、weather/，不要求重新传输约 7.64 GB 数据。

状态与不变项
------------
依据用户提供的旧版 ECS 记录，Ubuntu 24.04 / EPYC 9T95 / R 4.5.3 已跑过真实模型，
保留 14 个完整 Tile 检查点，其中 12 个 144 格点 Tile 有完整计时。
平均 350.88 s/tile、2.437 worker-s/cell；worker 观测峰值 4.70/4.84 GiB、主进程 1.75 GiB。
稀疏 pilot 的单 worker 约 0.5 GiB。以上是用户提供的观测，不是严格峰值上界。

本轮修改执行调度、监测和测试；原 12×12 Tile、气象包/pack_index、M0/M1、2000 个配对、
day_chunk=31、C++ 内核、支持域、热边界/逐日 CDD、科学检查点格式与 Tile-ID 汇总顺序保持不变。
新版尚未在目标 32/128 或 64/256 ECS 实测，本地合成测试不能证明真实服务器性能或绝对无 OOM。

环境与入口
----------
需要可恢复现有 brms 对象的 R 环境、C++ 工具链及：
ncdf4, digest, Rcpp, brms, posterior, rstan, dlnm, FNN, splines。
code/model_environment_versions.csv 是历史拟合记录，不是安装锁文件。
原始模型数值 QC 及每个 worker 的预测一致性检查仍执行。

config.R 保存现有科学配置；config_memory.R 继承它并覆盖执行参数与输出根目录。
run_server.sh 默认使用 config_memory.R。pilot 使用 config_pilot.R，输出 results_memory_pilot。
正式新版本输出 results_memory_v2，保留旧 results 及其 14 个基准检查点。

现有 ECS bundle 目录为 /data/era5/ERA5_workflow/03_server_bundle。
run_server.sh 自动优先通过 /data/era5/micromamba 的 era5 环境运行，无需交互式激活；
默认查找 $HOME/.local/bin/micromamba 或 PATH 中的 micromamba，兼容 nohup。
可用 MAMBA_ROOT_PREFIX、MAMBA_EXE、ERA5_MAMBA_ENV 指定已有安装；不自动安装任何环境。
其他机器已有合适的 R 时，可显式 ERA5_R_ENV=current 使用当前 PATH。
在 bundle 目录执行：
  cd /data/era5/ERA5_workflow/03_server_bundle
  bash run_server.sh env
  bash run_server.sh check
  bash run_server.sh full
  bash run_server.sh verify
可选小规模端到端检查（仍不是性能基准）：
  bash run_server.sh pilot

check 执行静态 SHA256、50 项数值回归、30 项 PSOCK/合成执行与恢复测试、
34 项内存调度测试（Linux 另加真实 /proc/cgroup 检查）、7 项启动入口测试、气象包验证和依赖/情景/资源预检。
full 自动进行真实模型和逐 worker 数值 QC，通过后才派发科学任务；不需要另开多轮人工性能试验。
如使用 tee 管道，先 set -o pipefail，保留失败退出码。

少量执行配置（config_memory.R）
------------------------------
scheduling="memory_lpt"：启用内存感知 LPT + Backfilling。
workers=64：初始进程池上限；CPU affinity/quota、RAM、剩余任务及初始化后内存可进一步限制。
autotune=FALSE：不再为比较 worker 数而分批跑基准任务。
memory$reserve_fraction=0.15：为系统及其他进程留至少 max(4 GiB,15%有效 RAM)。
memory$safety_factor=1.30：对 worker 峰值估计及增长预留使用安全系数。
memory$backfill=TRUE：允许较小任务利用内存空隙。
memory$max_rebuilds=2：全部任务排空、GC 后仍无法准入时，最多自动重建两次空闲进程池；0 可禁用。
progress_seconds=30：进度更新周期，memory_lpt 在无完成事件时也有轮询心跳。
resume="auto"：仅续接同版本且指纹一致的未完成运行。

32 vCPU/128 GiB 可将 workers 设为 32；64 vCPU/256 GiB 设为 64。
其余初始默认值相同。workers 是池上限，不是大型 Tile 的同时执行数量。
scheduling="lpt" 或 "fifo" 可用于对照；这两个模式没有动态内存准入，需另行保守设置 workers。

内存模型、准入与公平性
----------------------
初始完整 144 格点、31 天块、2000 draws、10 年、1 个 CDD 情景的未加安全系数估计为：
0.75 GiB 常驻 + 0.25 GiB 年度数组 + 4 GiB 瞬时矩阵 = 5 GiB；乘 1.30 后为 6.5 GiB。
年度部分随格点数、draws、年份和情景数变化，瞬时部分随格点数、day_chunk 和 draws 变化。
4 GiB 是由有限实测构造的保守约 60 个 double 矩阵等效额度，不是精确对象内存清单。
实际支持记录数写入计时/监测；平均支持率不能限制最密集时间块，所以内存仍按全支持估计。
观察到超出估计的独立增长只会上调内存模型；低观测不会自动下调。高驻留前序 Tile 不用于拟合小 Tile 的峰值。
耗时按小/中/大三组 seconds/cell 做简单平滑更新，不拟合复杂函数。

Linux 每秒左右读取当前 RSS、生命周期 VmHWM、MemAvailable、CPU affinity 与 cgroup v2 限制。
有效 RAM/可用量还受当前及祖先 memory.max、memory.high、memory.current 限制；CPU 遵守 cpu.max。
共享页面使 RSS 求和偏保守，VmHWM 属于进程生命周期，并非每个 Tile 的独立峰值。
仅当 HWM 比派发前增加时，才把这个新增上界辅助记入当前任务；采样峰值仍可能漏掉短暂尖峰。

每次派发都同时检查：
  预计主进程额度 + 所有 worker 承诺额度 <= 有效 RAM - 系统预留；
  主进程未来增长 + 所有 worker 尚未达到的峰值增长 <= MemAvailable/cgroup 余量 - 系统预留。
主进程额度至少 max(2 GiB,当前主进程 RSS × 安全系数)。
空闲 worker 的驻留内存仍计入；运行任务的预计峰值不因当前 RSS 较低而取消。
优先复用已有高驻留 worker，避免无谓扩张更多独立 R 堆。

主队列按预计耗时选任务；大 Tile 暂时不能准入时尝试小 Tile。
被插队达到 8 次或等待 60 秒的大 Tile 获得预留，固定其预计最早可启动时刻。
预留期间仍可回填：小任务预计耗时 × 1.5 + 1 秒必须落在该时刻之前，且把小任务峰值
作为可能保留的内存堆计入后，大任务在该时刻仍须有 worker 和足够内存。
预计释放额度基于当前保留内存和在途任务预测，不能保证真实完成时刻；每次实际派发仍执行实时内存准入。
预留时刻不会随等待而后移；超时或无法估计可用窗口时停止进一步回填，防止无限插队。
已有任务运行但没有新任务可准入时继续轮询和收取结果，保留连续队列，没有固定批次屏障。
完全空闲且仍无法准入时，仅向空闲 worker 请求一次 gc() 并重新测量。
若仍不满足预算，仅在全部任务已返回后关闭空闲池、确认进程退出，再以原池一半的上限重建。
单 worker 时可重建一个新进程。最多两次，每个新 worker 重新加载和执行预测 QC、登记 PID。
已提交 Tile 不重算；新池只接收本次尚未完成的 Tile。次数用尽、初始化失败或旧进程 10 秒内未退出则明确失败。
没有关闭正在计算的 worker，也不假定小 Tile 或 gc() 会使 R 释放前一块的大堆。

失败、中断与恢复
----------------
worker 初始化按顺序进行，每个初始化前复查预算；失败时清理已创建连接。
worker 计算错误、连接断开、回调/日志错误或监测不可用时，停止新派发，等待其他在途任务返回。
上述错误以及 STOP/中断均不触发自动重建；只恢复明确的全空闲内存准入失败。
监测失败不会静默降级到无内存约束的队列；有意使用 LPT 需要显式配置。
进程级强制退出仍依赖既有原子提交机制，不能保证未提交 Tile 被保存。

温和停止可在当前 run_*/STOP 创建空文件；调度器停止派发并排空。移除 STOP 后同命令恢复。
详见 DEPLOY_MEMORY.txt 的准确命令。已完成的运行不会自动重做。
检查点已经原子提交而 CSV/计时 sidecar 缺失时，恢复仍验证并跳过该 Tile，记录日志缺失即可。

03_run_ecs.R 继续保留执行层内容 MD5，并新增 11_memory_resources.R 的 MD5。
旧 14 个检查点不迁移、不修改指纹；本版本默认新结果根目录。
同版本中 workers、scheduling、memory、progress_seconds 不进入科学身份，允许中断后调整；
代码、数据、科学设置、后验配对、后端及 R/相关包环境变化仍接受严格指纹检查。

输出和审计
----------
results_memory_v2/run_*/global/：
  checkpoints/tile_XXXXX.rds：原有逐格年度配对后验科学检查点。
  progress.csv：有效完成计数、速率、ETA 估计依据和不确定性。
  dispatch_<session>.csv：候选队列排序。真正派发次序/完成以 memory_dispatch 为准。
  memory_dispatch_<session>.csv：派发/完成、Tile/格点、worker/PID、时间、状态、检查点大小；
    predicted_* 为预测；sampled_task_peak_rss 为采样/新增 HWM 辅助观测；
    process_lifetime_hwm_* 明确表示进程生命周期峰值，不冒充 Tile 独立峰值；
    committed/reserved_future/headroom 为派发时预算；reason 为 LPT/回填/预留原因。
  memory_monitor_<session>.csv：活跃/等待数量、被预留 Tile、可用内存、RSS 总和及未来增长预留。
    telemetry_valid 和 sample_time 区分当前有效采样与监测故障后的旧读数。
  resource_plan_<session>.rds：有效资源、初始池上限、保守参数和观测来源。
  memory_recovery_<session>.csv：每次成功重建的原因、前后池大小和剩余 Tile。
  worker_QC_<session>_rebuild_<n>.rds：重建池的逐 worker QC 记录；失败原因写入主日志。
  resume_<session>.csv：本会话已验证的旧检查点、格点数及计时 sidecar 是否存在。
  tile_timings_<session>.csv、tile_timings/tile_XXXXX.rds：原有分阶段时间与持久计时。
ETA 用剩余 Tile 大小分组成本和按时间加权的观测并行度，初期有明确不确定性；不扣除活动 Tile 已完成比例，
不包括最终汇总和未来外部内存压力变化。它是规划估计，不是运行时保证。

文件与科学配置
--------------
weather/packed_2010.nc ... packed_2019.nc、pack_index.rds、年度完成记录和格点清单均复用。
SHA256SUMS.txt 校验静态代码/说明/模型/pack_index；WEATHER_MD5SUMS.txt 校验年度 NetCDF。
config.R、config_memory.R、config_pilot.R 可编辑，运行时科学配置仍记录到指纹。

原 CDD 情景保留 theta0=18、slope=1/0.9683、kappa=1。
每日 CDD_M0=max(Ta_outdoor-18,0)、CDD_M1=max(Ta_outdoor-18-0.9683×ΔU,0)，先逐日截断再累计；
kappa=1 是室内变化等量转移到室外门槛的情景，截距不进入差值。
CDD 单位为 °C·day，面积加总为 km²·°C·day，不是耗电量。

下载整个 results_memory_v2 到本地 04_results，保留 global/checkpoints、grid_partitions、后验及运行元数据。
现有预处理确认开关被入口硬编码覆盖的问题在 DEPLOY_MEMORY.txt 单独记录，本次未改动这一历史科学行为。

监测/API 依据：
https://docs.kernel.org/filesystems/proc.html
https://docs.kernel.org/admin-guide/cgroup-v2.html
https://stat.ethz.ch/R-manual/R-devel/library/base/html/socketSelect.html
https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html
