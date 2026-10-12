# Thermal History Global — ERA5 workflow

使用已拟合的 M0/M1 模型，对全球 ERA5 数据计算热感觉边界变化及预先指定情景下的制冷度日（CDD），保留 2,000 个固定配对后验抽样。CDD 是气候指标，不直接等于耗电量。

范围为 2010—2019 年、原始 0.25° 网格、60°S—80°N；网格与陆地相交且 2015 年人口密度严格大于 1 人/km²。当前气象包包含 104,344 个格点、1,990 个 12×12 tile 和 3,652 天。

截至 2026-10-12，本地气象预处理和上传包验证已完成；真实模型的全域预测、目标 ECS 性能测试尚未完成。详细操作见 [本地说明](00_README.txt) 和 [服务器说明](03_server_bundle/README_SERVER.txt)。

## Git 保存什么

这个仓库保存源代码、科学配置、测试和必要的溯源说明。完整服务器上传包约 7.64 GB，需要在 Git 之外传输。克隆代码仓库不会取得模型、气象数据或现成的完整运行环境。

| 文件或目录 | Git 范围 | 原因 |
| --- | --- | --- |
| `.gitignore`、`.gitattributes`、本说明、`00_README.txt`、`RUN_LOCAL.ps1` | 保存 | 发布规则、操作说明和本地入口 |
| `03_server_bundle/code/*.R` | 保存 | 活动计算代码、原始计算定义和合成测试 |
| `03_server_bundle/code/model_environment_versions.csv` | 保存 | 历史模型环境记录；不是自动安装锁文件 |
| `03_server_bundle/config.R`、`config_pilot.R` | 保存 | 可审查的科学配置及试跑设置 |
| `03_server_bundle/run_server.sh`、`README_SERVER.txt` | 保存 | 服务器入口和运行说明 |
| `03_server_bundle/SHA256SUMS.txt`、`WEATHER_MD5SUMS.txt` | 保存 | 校验另行传输的代码、模型、索引和气象文件 |
| `02_local/mask/` 中的三个 R/PowerShell 脚本和 `README.txt` | 保存 | mask 的生成方法和核验方法 |
| `02_local/audit/check_month_reader.R`、`path_migration.csv` | 保存 | 数值核验脚本及现有 mask 溯源所需的迁移映射 |
| `90_archive/scripts/0.5°_ERA5-Figure1-population.R` | 保存，唯一归档例外 | 当前 mask 脚本仍从中提取人口面积分配函数 |
| `04_results/README.txt` | 保存 | 结果取回与保存要求 |
| `01_source_data/` | 忽略 | 原始气温、热历史和人口数据，约 50.49 GB |
| `02_local/` 其余内容 | 忽略 | mask 数据、海岸线、图表、审计产物、日志、R 库和临时目录 |
| `03_server_bundle/models/`、`weather/` | 忽略 | 模型及支持对象、NetCDF、pack_index 和格点/tile 数据表 |
| `03_server_bundle/results/`、`results_pilot/`、`logs/`、`UPLOAD_READY.txt` | 忽略 | 运行输出和只对本机已组装上传包成立的状态记录 |
| `04_results/` 其余内容、`90_archive/` 其余内容 | 忽略 | 大型结果、检查点、历史副本和拟合材料 |
| 凭据、R 会话和编辑器缓存 | 忽略 | 与代码版本无关的本机状态和私密配置 |

`90_archive` 中的人口脚本虽是历史文件，却仍是 mask 生成的依赖，因此不能整目录排除。原文件位置及内容保留，避免改变已有数据指纹。

忽略规则采用路径白名单：新建顶层文件或目录后，需要主动更新 `.gitignore` 才会进入候选提交；`code/` 下新建的 `.R` 文件会正常显示。被忽略的文件仍完整保留在本地。不要用 `git add -f` 批量绕过这些规则。

## 从克隆恢复可运行目录

新克隆只包含上述代码和说明。Git 不保存空目录，需要时创建 `models/`、`weather/` 等输入目录；服务器入口会创建 `logs/` 和 `results/`。

服务器预测至少需要另行取得：

- `03_server_bundle/models/` 下的 `fit_01-50-lag3.rds`、`fit_01-M0-base.rds`、`data_01-50-lag3.rds`、`pred_01-50-lag3.rds`。
- `03_server_bundle/weather/` 下的十年 `packed_YYYY.nc`、`pack_index.rds`、年度完成记录及格点/tile 清单。建议恢复整个已验证的 `weather/` 目录。
- 能恢复已有 brms 拟合对象的 R 环境、C++ 工具链，以及服务器说明列出的依赖。

从原始数据重新打包还需要 `01_source_data/` 和已验证的本地 mask；重新生成 mask 还需要人口栅格与 Natural Earth 海岸线，见 `02_local/mask/README.txt`。恢复旧 mask 的溯源核验可能需要迁移表指向的历史原文件；这些原文件并非全部随 Git 分发。

当前 `config.R` 中的本地根路径及 PowerShell 入口中的 R 安装路径是这台 Windows 机器的设置。换机器后按实际位置调整；不要把凭据写入这些已跟踪的文件。历史参考脚本中的旧路径作为溯源信息保留。

运行前必须明确填写 `cdd$scenarios` 的 `name/theta0/slope/kappa`，并完成气象预处理元数据核验及 `preprocessing_acknowledged` 设置。当前不预设 CDD 情景；默认 32 workers 也不代表已完成云端调优。

合成测试入口为 `RUN_LOCAL.ps1 -Mode test`，需要相应的本地 R 环境；通过合成测试不等于真实模型或 ECS 已验证。服务器使用 `bash run_server.sh check` 核验完整输入后，再进行试跑。

## 文件校验与换行符

`.gitattributes` 将常规文本和 Linux 入口固定为 LF，避免 Windows 的自动换行转换影响服务器运行。

`warm02_reference.R`、历史人口参考脚本和 `model_environment_versions.csv` 按原始字节保存。后者当前使用 CRLF，并已列入静态 SHA256 清单，不能为统一格式而改写。配置文件不在静态清单中，科学设置由运行指纹记录。

静态 SHA256 清单包含未随 Git 分发的模型和 `pack_index.rds`，所以在仅有代码的新克隆中不能期待完整校验通过。恢复数据后再核验；修改清单覆盖的代码或说明后，也必须同步更新对应 SHA256 条目。提交 Git 本身不会自动重算校验码。

## 后续发布

先在本地查看候选内容，再提交、选择远端可见性并发布。GitHub Desktop 中应重新添加现有仓库路径 `E:\ERA5\ERA5_workflow`；原来的 `Thermal-History-Global` 子目录已并入项目根目录，初始提交历史保留。

```powershell
git status --short
git add --dry-run .
```

检查候选内容后，才执行实际暂存并查看提交差异：

```powershell
git add .
git diff --cached --stat
git diff --cached
```

本次规则不使用 Git LFS：现阶段模型和数据通过完整上传包及校验清单单独管理。结果检查点仍需单独备份，Git 仓库不替代数据备份。

这里的“保存”表示技术上适合进入代码仓库，并不替代对代码归属和数据授权的判断；尚未添加开源许可证，也未选择远端的 Public/Private 可见性。
