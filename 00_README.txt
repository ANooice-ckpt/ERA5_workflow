ERA5 工作目录
=============
整理日期：2026-10-14
位置：E:\ERA5\ERA5_workflow
状态：本地预处理与全包验证已完成。十个 NetCDF 共 7.389 GB（6.881 GiB）；
原完整 03_server_bundle 上传目录约 7.643 GB（7.118 GiB）。用户提供的旧版 ECS 实测已完成 14 个真实 Tile。
本轮新增内存感知调度，已通过本地合成验证；新版尚未在 Linux/目标 ECS 实测，全域尚未完成。

先看这些
--------
02_local/mask/mask_global.png       全球格点分布图
02_local/mask/mask_summary.json     mask 数量及方法摘要
02_local/audit/LOCAL_PREPARATION_STATUS.txt  本地预处理最终状态
02_local/audit/weather_pack_by_year.csv    年度日期数、大小和校验码
02_local/logs/pack_optimized_2010_2019.log   打包续跑及最终验证日志
02_local/audit/PACKING_NOTES.txt     本地 I/O 优化、初始日志与切换记录说明
02_local/audit/PARALLEL_UPDATE_20261012.txt  连续队列、LPT 和计时更新及测试说明

文件分类
--------
01_source_data/     原始气温、热历史、人口数据。仅本地保留。
02_local/           mask、人口匹配核验、审计表、日志及本地 R 依赖。
03_server_bundle/   唯一的活动代码和服务器输入；上传时只传此目录。
04_results/         预留，用于接收服务器的正式计算结果。
90_archive/         旧代码版本、历史审阅、模型拟合全套材料；不参与日常运行。
RUN_LOCAL.ps1      本地打包/验收/回归测试入口。
BUILD_CODE_UPDATE.ps1  仅代码增量打包，不含模型/气象数据或服务器科学配置。

03_server_bundle 内部
--------------------
code/              流式计算和核验脚本，warm02_reference.R 为未改动的原始计算定义。
config.R           唯一主配置，供本地和服务器共同使用。
config_memory.R    继承科学配置，启用内存准入、回填和独立 results_memory_v2。
config_pilot.R     继承新执行配置，限制 300 个格点，输出 results_memory_pilot。
weather/           十年选定格点的 NetCDF、pack_index.rds、年度完成记录和格点/tile 清单。
models/            服务器实际需要的四个 RDS，保留原名。
results/           正式服务器结果预留目录。
results_pilot/     小规模试跑结果预留目录。
logs/              服务器日志预留目录。
run_server.sh      服务器入口；自动识别上传后的目录位置。
README_SERVER.txt  服务器操作说明。
DEPLOY_MEMORY.txt  Windows SCP 增量部署、正式运行、查看进度、中断和恢复的准确命令。

固定计算范围
------------
原始 0.25° 网格中心，南纬 60° 至北纬 80°。
网格与陆地/小岛有交集即纳入陆地候选；2015 年人口密度严格 >1 人/km²。
保留 104344 个格点；12×12 tile 共 1990 个非空任务。
2010—2019 年，共 3652 个 UTC 日；五项气象量保持源值，float64 存储。
2000 个后验配对位置保持不变；本地打包不加载模型，不进行全球预测。

本地重跑或核验
--------------
在本目录打开 PowerShell：
  .\RUN_LOCAL.ps1 -Mode verify
  .\RUN_LOCAL.ps1 -Mode test
  .\RUN_LOCAL.ps1 -Mode pack
pack 会核对原始输入与已完成年度文件；通过校验的年份可以跳过。
运行中不要更改原始数据、mask、打包设置或移动目录。

上传与取回
----------
首次部署上传整个 03_server_bundle；现有 ECS 更新只需 BUILD_CODE_UPDATE.ps1 生成的代码包，见 DEPLOY_MEMORY.txt。
服务器无需固定为 /data/era5，使用 bash run_server.sh 时自动设置 ERA5_BUNDLE_ROOT。
完成后，将服务器 results_memory_v2/ 的完整内容下载至本地 04_results/；旧 results/ 保留。
保留 run_*/global/checkpoints 的逐格年度后验、pack_index 和运行元数据，
以后才能按区域重新聚合，不能只保存 PNG 或分位数表。

CDD 情景与运行前事项
--------
config.R 已按 2026-10-12 的确认设置 theta0=18°C、slope=1/0.9683、kappa=1。
图片给出 ΔTa≈0.9683×ΔU；现有公式使用 kappa/slope，因此 slope 填倒数。
kappa=1 是将室内温度变化等量转移到室外 CDD 门槛的情景假设；截距 0.583 在差值中抵消。
每日 CDD_M1=max(Ta_outdoor-18-0.9683×ΔU,0)，逐日截断后再累加。
仍需完成预处理元数据核验及真实模型数值 QC；历史确认开关被入口覆盖的问题单独记录，未混入本次调度改造。
新默认 config_memory.R：workers=64 是池上限；按 CPU、RAM 和实测 RSS 自动限制实际派发。
调度为 Memory-Aware LPT + Backfilling；留至少 15%/4 GiB 系统余量，峰值安全系数 1.30。
大 Tile 被插队 8 次或等待 60 秒后预留资源；所有在途任务未来内存增长、空闲 worker 驻留均计入预算。
服务器运行后，global/tile_timings_<session>.csv 可查看各阶段耗时，
global/dispatch_<session>.csv 记录待派发顺序，progress.csv 记录完成进度。
memory_dispatch_<session>.csv 记录实际准入/完成，memory_monitor_<session>.csv 记录内存及队列状态。
本地 test 包含 50 项数值回归、27 项实际 PSOCK/合成执行恢复、28 项内存调度检查。
Linux 真实监测由服务器 check 自动检查；本机没有 Linux/WSL，不将合成结果当成新版 ECS 实测。

迁移与溯源
----------
路径迁移清单：02_local/audit/path_migration.csv。
原始数据在 E: 同盘迁移，没有重新生成或替换。
历史 RDS 内记录的旧路径作为溯源记录保留；mask 验证器通过迁移表解析新位置。
模型拟合全套材料留在 90_archive/model_fitting，上传包只包含所需四个模型/支持文件副本。
v1.1 原发行包及 SHA256 留在 90_archive；活动配置已适配当前目录。
