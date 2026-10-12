独立空间 mask 测试（不运行原脚本主程序，不计算模型预测）
网格中心纬度范围：-60 至 80；网格间隔：0.25°；人口年份：2015。
全部格点 807,840；陆地格点 260,201；确认 >1 人/km² 的陆地格点 104,344。
确认 ≤1 的陆地格点 155,857；人口数据不足 0。
保留陆地格点的 40.1013%；格点数量减少 59.8987%（不代表运行时间同比减少）。
方法：读取新版 warm02 的 config 和原始坐标。原中心陆地规则仅作对照；主 mask 改用 Natural Earth 1:10m 陆地及小岛多边形，与网格有交集即入选陆地候选。
仅复用人口脚本中可接受任意目标分辨率的保守面积分配函数，不运行其半度主程序。
人口原始栅格值是每个像元的人数，不是人口密度。依据真实像元边界和 WGS84 椭球面积分配人数，除以完整 0.25° 网格面积。
没有按网格内陆地面积校正分母。海岸网格也使用完整网格面积。
陆地使用 terra::rasterize(touches=TRUE)，跨日期变更线网格合并两半；不是 ERA5 官方陆海比例 mask。结果受边界数据分辨率影响，并非保证包含世界上每一个微小岛礁。
绘图用真实矢量海岸线遮去海面，保持原始 0.25° 计算网格；没有对科学数据插值或平滑。
NoData 和源数据范围外区域不视为零人口。已知人口贡献足以令密度 >1 的格点可确认保留；否则人口覆盖不全时标为未知。
严格使用 >1，不是旧人口脚本的 >=1。气候/模型逐日支持、有效天数等条件尚未检查。
NetCDF 保留源经纬度排序。mask_status: 0 无陆地交集；1 保留；2 低密度排除；3 人口未知。
selected: 1 保留；0 排除/海洋；缺失为人口未知。RDS 和 CSV 保留原生产 cell_id。
Windows rerun: run_mask.ps1 (sets UTF-8 locale and temporary directory).
Windows 若 R 无法创建临时目录，请先将 TMPDIR、TMP、TEMP 设为本测试目录下的 _tmp。
依赖 ncdf4、maps、terra、sf、ggplot2、ragg、jsonlite；本次新增包安装在 ../runtime/Rlibrary，不修改系统库。
人口数据说明：https://human-settlement.emergency.copernicus.eu/documents/GHSL_Data_Package_2023.pdf
海岸线来源：https://www.naturalearthdata.com/downloads/10m-physical-vectors/
