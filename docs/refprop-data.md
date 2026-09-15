# REFPROP 离线建表与数据边界

## 使用方式

普通计算只需将工程根目录加入 MATLAB 路径。`+gasprop` 读取 `data/` 中的 MAT 文件，不加载 REFPROP，也不添加或运行 `tools/`、`validation/`。

只有重建数据时显式执行：

```matlab
addpath('tools');
outputs = buildRefpropTables;
% 可指定安装位置和网格密度：
% outputs = buildRefpropTables(InstallPath="D:\REFPROP", ...
%     TemperatureCount=241, VolumeCount=641, PressureCount=241);
```

生成结果：

- `data/{helium,hydrogen,nitrogen,argon,air,carbon-dioxide}-refprop.mat`：变量 `tabularGas`，六种工质通过 `gasprop.table.build` 编译的运行数据。
- `data/rk-reference.mat`：变量 `rkReference`，He、H2、N2、Ar、Air、CO2 的一维参考物性。
- 默认不保存原始采样或结果日志；显式 `SaveSource=true` 或 `Compile=false` 才输出 `*-refprop-source.mat`。

新增工质的 REFPROP 配置入口是 `tools/+refpropbuild/fluidConfig.m`；核心定义同时加入 `gasprop.gasCatalog` 和 `GasProperties`。建表函数会检查其 REFPROP 模型和熔化线；未提供可用相界时终止，不猜测固相边界。通用 `gasTable` 为全部六种工质生成二维 EOS 和输运表。可用 `buildRefpropTables(Gases=["CO2","Air"])` 只重建选定二维表，一维参考仍保持全部工质。

## 本机源数据与接口

实际使用 `D:\REFPROP\REFPRP64.DLL`，DLL 查询版本 **9.1**，氦气源文件 `helium.fld` 使用 Ortiz-Vega 等的 2013 年中间版 EOS。它与 REFPROP 10 的后续氦气方程不能视为完全相同。

用户提供的 `refpropm.m`、`rp_proto.m`、`rp_proto64.m` 原样保存在 `validation/external/refprop/`。本机缺少旧接口所需的 thunk DLL；离线初始化器使用 [NIST 官方 9.1.1 thunk 文件](https://trc.nist.gov/refprop/FAQ/MATLAB/9.1.1/REFPRP64_thunk_pcwin64.dll)，存于 `tools/+refpropbuild/`。它仅将旧原型的 thunk 符号名称映射到 `int32/cstring` 名称，保留用户指针参数类型。没有修改 REFPROP 安装或用户接口。

该补充文件是 [NIST 官方 MATLAB 接口说明](https://github.com/usnistgov/REFPROP-wrappers/blob/master/wrappers/MATLAB/legacy/README.rst) 列出的依赖。其 SHA256 为 `0f1ab54186e3287fb064184f374b15b0c10ac1a350a80d97907223502d491530`。实际源 DLL、各流体文件、用户接口和 thunk 的 SHA256、生成时间、MATLAB 版本都写入生成文件的 `SourceMetadata`。

## 二维表

氦气表的运行范围为 **10–800 K、0.1–20 MPa**。REFPROP 实测氦气临界温度 5.1953 K，20 MPa 下熔点 5.152163725 K，所需区域均位于气体或超临界单相区；其余五种工质的温区见下方表格。

EOS 使用对数温度轴和对数比容轴。默认 161 个温度点另加入支持温区内的参考温度；氦气共 162×481 节点。其他工质额外使用对数 `(T-Tc)` 温度节点和临界密度附近的比容节点，加密单相临界附近的陡峭变化；准确尺寸见 MAT 文件的 `SourceMetadata.EOSGridSize`。所有拼接网格先对齐精确端点，再按 64 倍机器精度合并近重复点，避免 `800` 与 `799.9999999999998` 形成极薄插值单元并破坏导数。比容覆盖整个 T-p 请求包络并添加 0.01% 余量。采样包含：

- `Z=p/(rho*R*T)`、内能、熵、Cv；
- `dp/dT|rho` 和 `dp/drho|T` 的 REFPROP 导数；
- 零密度理想比热参考值；
- 完整来源、节点数、运行边界和每个温度的最高内部压力及该压力下熔点。

`R` 使用运行时所选气体常数；因此不直接保存 REFPROP 自身定义的 Z。`refpropm` 压力及 `#` 导数由 kPa 转为 Pa；其 `R` 输出为 `d(rho)/dP`，取倒数后再乘 1000 得到 `dp/drho`。

T-v 矩形高温高密度角点超过实际请求压力，仅作为插值内部覆盖。He/H2/N2/Ar/Air/CO2 的最高内部压力约为 **705/540/439/450/442/313 MPa**，均低于各自 EOS 上限 **1000/2000/2200/1000/2000/800 MPa**。公开工作域最高仍为 20 MPa；温压范围保护始终开启，不受可选校验开关影响。

每个温度先求最大密度节点的压力，再以 `MELTP(p_max(T))` 求该压力下的熔点，要求熔点在流体文件 `#MLT` 声明范围内且低于当前温度。此方法避免在 H2 熔化线 400 K、Ar 熔化线 700 K 声明上限以外强行调用 `MELTT(800 K)`。同时检查 EOS 最大密度、压力、每节点单相质量分数及稳定性量。

本机氦气输运模型上限为 **100 MPa**，不能给完整 T-v 矩形强行填满黏度与导热数据。故全部工质的 `TransportTable` 单独使用 **T-p 网格**，只采样 0.1–20 MPa 合法请求区，并检查各自输运压力及密度上限。氦气默认 162×161 节点；其余工质还在临界压力附近加密。输运使用双线性插值，为解析近临界导热峰，H2/N2/Ar/Air 的温度与对数压力单元各二等分，CO2 各四等分；因子由 `fluidConfig.TransportRefinement` 配置。精确尺寸见 `SourceMetadata.TransportGridSize`。核心根据已经求出的压力查询此表，不需要 REFPROP。

模型上下限来自本机 9.1 的 `LIMITSdll('EOS'/'ETA'/'TCX')`：氦气 EOS 为 2.1768–2000 K、1000 MPa；黏度与导热均为 2.1768–1500 K、100 MPa。`LIMITS` 只返回模型范围，相态还由 `MELTP`/`MELTT` 熔化线及临界温度单独检查。氦气文件的熔化线为 McCarty–Arp (1990) 模型，声明温度范围 2.1768–1500 K；所用 10–800 K 均在其范围内。这里的范围以安装文件为准，未套用其他 REFPROP 版本的参数。

`buildRefpropTables` 显式临时开启 `gasprop.validation(true)`，正常结束或报错后恢复原开关和 path。范围与相态检查统一位于 `validation/+gaspropcheck/refpropBuild.m`；整表一致性审计和独立误差检验仍需在 `validation/` 显式运行。

## RK 参考物性和范围

`THERM0dll` 计算方程中的理想气体项 `cp0(T)`，转换为质量基准。`TRNPRPdll(T,D=0)` 返回零密度 `mu0(T)`、`k0(T)`；没有用高压真实 Cp 代替理想 Cp。函数定义与单位见 [REFPROP Legacy API](https://refprop-docs.readthedocs.io/en/latest/DLL/legacy.html)。

| 工质 | REFPROP 文件 | 参考表温度范围 K | 物理 Tc K | 20 MPa 熔点 K |
|---|---|---:|---:|---:|
| He | helium.fld | 10–800 | 5.1953 | 5.15216 |
| H2，正常氢 | hydrogen.fld | 34–800 | 33.145 | 19.46710 |
| N2 | nitrogen.fld | 127–800 | 126.192 | 67.43720 |
| Ar | argon.fld | 151–800 | 150.687 | 88.71659 |
| Air，干空气伪纯流体 | air.ppf | 133–800 | 132.6312 | 63.24201 |
| CO2 | co2.fld | 305–800 | 304.1282 | 220.67705 |

上述温区也用于全部 Tabular 表，最高压力均为 20 MPa。参考表默认每个工质 401 个对数温度点，参考温度取 `max(300 K,Tmin)`：CO2 为 **305 K**，其余为 300 K。不将 300 K 强行加入 CO2 气态 EOS 或参考运行温区。温度下限保守地高于物理临界温度、20 MPa 熔点及模型温度下限；这保证整个 0.1–20 MPa 压力区间不会跨入亚临界液相、两相或固相，允许超临界单相流体。低压下本可使用的部分较低温气態区因此也被保守排除。温压范围保护始终开启，表外查询明确报错；`gasprop.validation(false)` 仅关闭额外诊断，不关闭范围保护。

CO2 使用本机 `M=0.0440098 kg/mol`、`Rmol=8.31451 J/(mol K)`、`R=188.92405782348476 J/(kg K)`。305 K 零密度值为 `Cp0=851.01298552020774 J/(kg K)`、`mu0=1.5255769724942452e-5 Pa s`、`k0=0.017152634525186769 W/(m K)`；300 K 理想极限值分别为 `845.84601038795631`、`1.5013923982513168e-5`、`0.016747259697437032`。后者可供 Ideal 常物性定义，但不表示 300 K 全压力域为气态。

Sage 拟合 RK 的有效临界点可能高于物理临界点，RK 运行时会进一步收紧有效下限，详见 `sage-rk-comparison.md`。参考表的温度范围并不覆盖被排除的低温相变区。

来源元数据同时记录 REFPROP 质量气体常数和本代码的运行时气体常数，便于审核旧配置取值与 REFPROP 的差异。`cp0` 是 REFPROP 原始质量基准，RK 热力学残差和理想部分则使用本代码同一 `R`；不得把这组数据宣称为 REFPROP 完整状态方程。

原 H2 配置 `R=4157.2 J/(kg K)` 比本机 REFPROP 的 `4124.487568704486` 高约 0.793%；已统一修正为 `M=0.00201588 kg/mol`、`R=8.314472/M`。He 保持 NISTHelium 使用的公共 R；N2、Ar、Air 保留现有小幅舍入值，并在来源元数据明确记录。
