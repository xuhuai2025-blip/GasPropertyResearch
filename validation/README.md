# 独立校验

本目录保存保护性检查、整表审计和 REFPROP 对照接口。没有自动运行的测试或结果文件。

## 会话开关

```matlab
addpath("D:\VibeCodingWorkSpace\GasPropertyResearch")
enabled = gasprop.validation(); % 新会话默认 false
gasprop.validation(true);       % 开启检查并加入 validation 路径
gasprop.validation(false);      % 关闭检查
```

`+gaspropcheck/` 统一存放输入、范围、数据结构和状态检查。开关作用于当前 MATLAB 会话中的全部对象；构造或加载数据的检查需要在这些操作之前开启。关闭时调用者负责提供有效数据和工作域内状态，不保证越界或无效状态被拒绝。求根收敛处理和模型分派错误不受此开关影响。

开关只控制保护性检查，不会启动下面的整表扫描或外部对照。

## 整表审计

```matlab
root = "D:\VibeCodingWorkSpace\GasPropertyResearch";
addpath(root,fullfile(root,"validation"))
gas = gasprop.GasProperties.fromName("He");
table = gasprop.table.load(fullfile(root,"data","helium-refprop.mat"),gas);
audit = validateTabularTable(table);
% 包含完整 EOS 插值支撑区：
audit = validateTabularTable(table,struct("IncludePadding",true));
```

入口直接评估原生 T-v 曲面，每个单元双向各取 3 个内部点，分块统计稳定性与热力学恒等式。它保持原开关和 path 不变，不借助生产状态范围检查筛掉填充点。默认只统计工作压力域；`IncludePadding=true` 包含完整插值支撑区。

主要验收条件：

- 压力、Cv、Cp、压力对密度导数和声速平方有限且为正。
- 能量恒等式残差 `abs(u_v-(T*p_T-p))/max(abs(T*p_T),abs(p)) <= 0.001`。
- 熵导数的缩放误差不超过 0.05。
- 原 u_v 相对指标保留为 `EnergyVolumeDerivative` / `LegacyRelativePassed`，便于检查两个压力项相消造成的相对误差放大。

入口只返回报告，不保存文件；`ErrorOnFailure=true` 可在审计失败时抛错。其他参数为 `Fractions`、`Tolerance`、`EnergyIdentityTolerance` 和 `TemperatureChunkSize`。

## REFPROP 对照

```matlab
report = runRefpropValidation();
report = runRefpropValidation(struct("Gases",["N2","CO2"],"IncludeRK",false));
% 只有显式指定输出目录才保存 MAT/JSON：
report = runRefpropValidation(struct("OutputFolder",fullfile(tempdir,"gasprop-audit")));
```

入口临时开启 `gasprop.validation(true)`，加入离线工具及用户的 `external/refprop/` 接口；正常结束或报错均恢复原开关和完整 path。默认安装位置为 `D:\REFPROP`，可用 `RefpropPath` 覆盖。

对照使用六种工质各自数据域内的独立 T-p 网格，包含边界与近临界点；内域同时按同 T、同 REFPROP 密度查询。Tabular 密度、Cv、Cp、声速、Z、压力导数和输运的默认最大相对误差门槛为 0.5%。u/h 另用热容量尺度验收，原始 u/h/s 及误差位置完整保留；Z 使用相同运行库 R，REFPROP 原生 Z 另存。

RK 只报告近似误差，不适用 Tabular 的精度门槛；同 REFPROP 密度下，RK 预测压力可能超出其范围，相关失败分别保留。采样对照不构成整个连续域的误差上界证明。

其他选项：`TableFile`（只选一种工质）、`ReportName`、`TemperatureK`、`PressurePa`、`RelativeTolerance`、`DeepTableAudit`。完整原始采样和对照结果应放在调用者选择的输出目录，默认不写入工作目录。
