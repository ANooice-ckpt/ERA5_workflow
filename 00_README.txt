ERA5 工作目录
=============
整理日期：2026-10-12
位置：E:\ERA5\ERA5_workflow
状态：本地预处理与全包验证已完成。十个 NetCDF 共 7.389 GB（6.881 GiB）；
完整 03_server_bundle 上传目录约 7.643 GB（7.118 GiB）。服务器模型预测尚未运行。

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

03_server_bundle 内部
--------------------
code/              流式计算和核验脚本，warm02_reference.R 为未改动的原始计算定义。
config.R           唯一主配置，供本地和服务器共同使用。
config_pilot.R     读取主配置，限制为 300 个格点并输出到独立 results_pilot。
weather/           十年选定格点的 NetCDF、pack_index.rds、年度完成记录和格点/tile 清单。
models/            服务器实际需要的四个 RDS，保留原名。
results/           正式服务器结果预留目录。
results_pilot/     小规模试跑结果预留目录。
logs/              服务器日志预留目录。
run_server.sh      服务器入口；自动识别上传后的目录位置。
README_SERVER.txt  服务器操作说明。

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
上传整个 03_server_bundle，保持其内部结构；不需要上传 01_source_data、02_local、90_archive。
服务器无需固定为 /data/era5，使用 bash run_server.sh 时自动设置 ERA5_BUNDLE_ROOT。
完成后，将服务器 results/ 的完整内容下载至本地 04_results/。
保留 run_*/global/checkpoints 的逐格年度后验、pack_index 和运行元数据，
以后才能按区域重新聚合，不能只保存 PNG 或分位数表。

尚待确定
--------
config.R 的 cdd$scenarios 仍为空；需要运行前确定 theta0、slope、kappa。
这不影响本地气象打包；服务器正式预测前必须填写。
workers 目前为 32，未依据真实 ECS 性能调优；本地输出不代表云端已跑通。
调度现为连续动态队列＋LPT；默认按 tile 选中格点数从多到少派发。
服务器运行后，global/tile_timings_<session>.csv 可查看各阶段耗时，
global/dispatch_<session>.csv 记录待派发顺序，progress.csv 记录完成进度。
本地 test 入口包含原有 48 项数值回归和 24 项实际 PSOCK 合成并行测试。

迁移与溯源
----------
路径迁移清单：02_local/audit/path_migration.csv。
原始数据在 E: 同盘迁移，没有重新生成或替换。
历史 RDS 内记录的旧路径作为溯源记录保留；mask 验证器通过迁移表解析新位置。
模型拟合全套材料留在 90_archive/model_fitting，上传包只包含所需四个模型/支持文件副本。
v1.1 原发行包及 SHA256 留在 90_archive；活动配置已适配当前目录。
