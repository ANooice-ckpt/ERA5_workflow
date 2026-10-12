ERA5 server bundle
==================
将整个目录上传到 Linux ECS 的数据盘，保持 code / models / weather 结构。
活动科学配置是 config.R。脚本自动使用当前 bundle 目录。

准备
----
需要可恢复现有 brms 拟合对象的 R 环境、C++ 工具链，以及以下 R 包：
ncdf4, digest, Rcpp, brms, posterior, rstan, dlnm, FNN, splines。
模型来源的包版本记录在 code/model_environment_versions.csv；它是历史环境记录，
不是自动安装锁文件。原始拟合和服务器预测的数值 QC 仍由 03_run_ecs.R 执行。

在 config.R 中预先填入 CDD 情景的 name/theta0/slope/kappa。
当前 workers=32；根据实际内存、单 worker 峰值和试跑速度调整。
默认 scheduling="lpt"：按 tile 中选中格点数降序派发，耗时相同时按 tile ID 排序。
全部待算 tile 共用连续动态队列，空闲 worker 立即领取下一块，没有 2×workers 的批次屏障。
这只是耗时估计；实际时间还取决于支持比例、CDD 情景数量和 I/O。
可设置 scheduling="fifo" 做对照；两者均保留原 tile ID、固定抽样位置和最终 tile-ID 汇总顺序。
LPT 会先同时运行大块，调整 workers 时应观察这段的内存峰值；每个 worker 的 BLAS/OpenMP 仍为 1 线程。

命令（在本目录执行）
--------------------
  bash run_server.sh check
  bash run_server.sh pilot 2>&1 | tee logs/pilot.log
  bash run_server.sh full 2>&1 | tee logs/full.log
  bash run_server.sh verify

check 检查静态文件 SHA256、合成数值回归、实际双 worker PSOCK 调度测试、十年气象包 MD5、模型文件存在性、
依赖和情景配置。pilot 使用相同科学设置和 2000 个抽样位置，只计算前 300 个格点，
结果写入 results_pilot；它用于端到端检查，不能当成全球结果或代表性速度估计。
full 写入 results。shell 管道若用于自动作业，请先执行 set -o pipefail，保留失败退出码。

续算
----
默认 resume="auto"：仅恢复未完成的同一版本运行；已完成运行不会重复启动。
改变情景、代码、mask、tile 或其他科学设置后，使用 resume="new" 开始新运行。
旧 v1 中途断点不能强行套用 v1.1 引擎。
本次并行更新也改变了代码指纹，不能强行续接更新前的运行；现有本地气象包可直接复用。
在本次版本内，仅改变 workers、scheduling 或 progress_seconds 不改变科学设置，可以续接同一未完成运行。
某个 worker 报错后停止派发新任务，等待已派发任务返回并保留成功检查点，然后报告错误。
进程崩溃或强制中断时，仍以已原子提交且验证通过的检查点为准。
config_pilot.R 继承 config.R，仅修改格点上限和输出目录。

运行进度与性能计时
------------------
results/run_*/global/progress.csv：已完成 tile/格点、吞吐量和剩余计算时间估计。
progress_seconds=30 控制完成事件触发的更新频率；长时间没有 tile 完成时不会凭空更新进度。
LPT 按大块优先运行，且不同地区支持比例不同，初期 ETA 只是估计；最终汇总时间不包含在 ETA 中。
dispatch_<session>.csv：本次会话待派发队列，含 tile ID、派发排名和格点数成本估计。
  遇到错误后，队列尾部可能尚未派发；实际完成情况以计时表和有效检查点为准。
tile_timings_<session>.csv：每完成一块追加一行，包含 worker/PID、格点数、预测行数、
  检查点字节数，以及读取、支持判断、预测、变换、CDD/warm 累加、汇总、检查点写入及回读校验的秒数。
tile_timings/tile_XXXXX.rds：每块成功提交后的持久计时记录，含运行指纹。
  分段时间是 worker 的墙钟时间，各段互不重叠；并行 worker 的秒数之和不是整个作业耗时。
  total_seconds 包含气象文件关闭与检查点提交，不含随后写入这份小型计时报表的时间。
  旧 cp$elapsed_seconds 为提交前计时，性能分析应使用新的计时报表。
checkpoints/tile_XXXXX_error.rds：失败块的阶段、错误和截至失败时的计时。
计时属于运行审计，不替代科学检查点，也不参与数值汇总。

本次实现只涉及连续调度、LPT 和分段计时。worker 初始化、C++ 矩阵融合、kNN 搜索树复用、
最终汇总 I/O 的进一步优化，留待目标 ECS 的分段计时确定瓶颈后处理。

文件与校验
----------
weather/packed_2010.nc ... packed_2019.nc：原始 0.25° 气象值的稀疏提取，不插值。
weather/pack_index.rds：格点索引、完整格点面积、人口、tile、日期和输入溯源。
weather/packed_YYYY.done.rds：年度完成记录；selected_cells.csv / tiles.csv 便于浏览。
SHA256SUMS.txt：代码、启动脚本、模型、pack_index 的内容校验。
WEATHER_MD5SUMS.txt：十年 NetCDF 的内容校验，可执行 md5sum -c WEATHER_MD5SUMS.txt。
可编辑的 config.R / config_pilot.R 不列入静态 SHA256，实际运行配置会记录在科学指纹内。

输出与取回
----------
results/run_*/global/：全域汇总、元数据、grid_partitions、checkpoints。
下载整个 results 到本地 ERA5_workflow/04_results；不要遗漏 checkpoints。
全球 CDD schema 为 v2：无合格数据的期次输出 NA 和 no_eligible_support 状态。
CDD 单位为 °C·day，面积加总单位为 km²·°C·day；不是电力消耗。
真实模型和 ECS 性能试跑尚未由此次本地预处理完成。
