# GasPropertyResearch

可独立使用的 MATLAB 气体物性包，提供统一的状态与输运接口。正常计算只需将项目根目录加入 path。

```matlab
root = "D:\VibeCodingWorkSpace\GasPropertyResearch";
addpath(root)

gas = gasprop.create("He","tabular");
s = gas.stateTP(77,8e6);          % T [K], p [Pa]
rho = gas.density(77,8e6);        % kg/m^3
[mu,k,pr] = gas.transport(77,8e6); % Pa s、W/(m K)、无量纲
```

## 工质与模型

内置 He、H2（正常氢）、N2、Ar、Air（固定组成伪纯流体）和 CO2。工质别名见 `gasprop.gasCatalog()`，模型名称见 `gasprop.modelCatalog()`。

| 模型 ID | 名称 | 算法与支持 |
|---|---|---|
| `ideal` | Ideal Gas | 常比热理想气体；六种工质 |
| `rk` | RK Gas | RK 方程及偏差函数；六种工质，附离线零压 Cp 和输运参考曲线 |
| `tabular` | Tabular Gas | 离线编译的 T-v 热力学曲面与 T-p 输运表；六种工质均附 REFPROP 数据 |
| `helmholtz` | NIST He-4 Helmholtz | 氦气专用基本方程；别名 `nist-he4-helmholtz` |

`GasProperties` 是 value class，修改模型应接收返回值：

```matlab
gas = gasprop.GasProperties.fromName("N2");
gas = gas.withEquationOfState("rk");
gas = gas.withTransportModel("reference");
```

输运可选 `sutherland`、`nist-helium`、`reference` 和 `tabular`。RK 默认使用离线零密度参考曲线；Tabular 默认使用表内输运数据。RK 的 `state.Cp` 是状态热力学比热，`idealCp(T)` 是零压参考比热，两者不可混用。[RK 与 SAGE 说明](docs/sage-rk-comparison.md)给出公式与系数。

## 接口与适用范围

- `gas.state(T,rho)`：温度 K、质量密度 kg/m³。
- `gas.density(T,p)` / `gas.stateTP(T,p)`：温度 K、压力 Pa。
- `gas.transport(T,p)` / `gas.transportAtDensity(T,rho)`：返回黏度、导热系数和 Pr。
- `gas.equationOfStateMetadata()`：模型与数据范围。

状态字段包含 Pressure、InternalEnergy、Enthalpy、Cv、Cp、PressureTemperatureDerivative、PressureDensityDerivative、SoundSpeed、Z、R、Density，均采用 SI 和质量基准。Entropy 为 RK/Tabular 的扩展字段；不同模型的绝对焓、熵参考态不保证一致。对象的 Cp/Cv 属性为工质常量，实际状态比热使用 `state.Cp/state.Cv`。

随包 Tabular 与 RK 数据均覆盖 **0.1–20 MPa**：

| 工质 | Tabular 温度 K | RK 温度 K |
|---|---|---|
| He | 10–800 | 10–800 |
| H2 | 34–800 | 约 35.29–800 |
| N2 | 127–800 | 127–800 |
| Ar | 151–800 | 151–800 |
| Air | 133–800 | 133–800 |
| CO2 | 305–800 | 305–800 |

精确范围以元数据为准。这里采用整个压力区间都适用的保守矩形域，会排除部分低温低压气态工况。NIST He EOS 的独立范围为 20–1500 K，旧输运相关式为 20–830 K、密度不超过 160 kg/m³；10 K 氦气可使用 Tabular。

## 校验开关：默认关闭

所有保护性输入、范围、数据结构和状态检查统一放在 `validation/+gaspropcheck/`，由当前 MATLAB 会话的一个开关控制：

```matlab
gasprop.validation()       % 查询开关；新会话默认 false
gasprop.validation(true)   % 开启，并将 validation 加入 path
gasprop.validation(false)  % 关闭
```

开关影响本会话的全部气体对象。需要检查构造或表加载时，应先开启再创建对象。关闭时由调用者保证输入、数据和温压域有效，**不保证越界拒绝或无效状态诊断**；必要的数值求根收敛处理和模型分派错误仍保留。

开启开关不会自动运行整表扫描或 REFPROP 对照。需要独立校验时显式调用：

```matlab
addpath(fullfile(root,"validation"))
table = gasprop.table.load(fullfile(root,"data","helium-refprop.mat"), ...
    gasprop.GasProperties.fromName("He"));
audit = validateTabularTable(table); % 不需要 REFPROP
comparison = runRefpropValidation(); % 需要本机 REFPROP
```

REFPROP 对照会临时开启保护性检查，正常结束或报错后恢复原开关和 path；整表审计直接评估曲面，保持开关不变。两者默认只返回内存报告，显式设置 `OutputFolder` 才保存 REFPROP 结果。[校验说明](validation/README.md)列出参数。

## 目录与建表

```text
+gasprop/
  create.m, GasProperties.m, gasCatalog.m, modelCatalog.m, validation.m
  +eos/                 四种 EOS 与内部数值算法
  +transport/           输运模型
  +table/               表加载与编译
data/                   六种运行表及 RK 参考曲线
tools/                  可选 REFPROP 离线建表工具
validation/
  +gaspropcheck/        由开关控制的保护性检查
  validateTabularTable.m
  runRefpropValidation.m
  external/refprop/     用户提供的外部接口
docs/                   方程与数据来源说明
```

核心求值不调用 REFPROP 或离线工具。加载已编译 MAT 表会复用已有曲面；自定义原始表用 `gasprop.table.build(raw,gas)` 编译，也可向 `gasprop.create("N2","tabular",source)` 提供 struct、MAT 或 SAGE 文本表。

重新生成随包数据：

```matlab
addpath(fullfile(root,"tools"))
buildRefpropTables();
% buildRefpropTables(Gases=["CO2","Air"]);
```

此步骤需要 REFPROP。`buildRefpropTables` 会显式临时开启校验，正常结束或报错后恢复原开关和 path；范围与相态检查位于 `validation/+gaspropcheck/refpropBuild.m`。默认保存运行表与 RK 参考曲线；`SaveSource=true` 才保留原始采样数据。详见 [REFPROP 数据说明](docs/refprop-data.md)。
