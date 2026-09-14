# GasPropertyResearch

可独立复用的 MATLAB 气体物性包。它从当前 SIMPLE 研究工程中提取，但不依赖求解器、项目 JSON、GUI 或多控制体闭合逻辑。将本目录加入 MATLAB path 后即可使用。

```matlab
root = "D:\\VibeCodingWorkSpace\\GasPropertyResearch";
addpath(root)

he = gasprop.create("He");
s = he.stateTP(77,8e6);                 % T [K], p [Pa]
[mu,k,pr] = he.transport(77,8e6);       % Pa s, W/(m K), -

n2rk = gasprop.create("N2","RK Gas");
srk = n2rk.state(300,5);                % T [K], rho [kg/m^3]
```

`GasProperties` 是 value class。选择 EOS 或输运模型时，应接收返回值：

```matlab
gas = gasprop.GasProperties.fromName("N2");
gas = gas.withEquationOfState("tabular",myTable);
gas = gas.withTransportModel("sutherland");
```

## 可选气体模型

`gasprop.modelCatalog()` 返回下表对应的结构化目录。模型 ID、展示名均可传给 `withEquationOfState`；`helmholtz` 与 `nist-he4-helmholtz` 是同一 He 专用模型的 ID/别名。

| ID | 名称 | 核心算法 | 适用范围 |
| --- | --- | --- | --- |
| `ideal` | Ideal Gas | `p = rho R T`，常比热 | He、H2、N2、Ar、Air |
| `rk` | RK Gas | Redlich-Kwong 立方状态方程及一致的偏差函数 | 纯 He、H2、N2、Ar；仅超临界温度、单相稳定区域 |
| `tabular` | Tabular Gas | SAGE `Z(v,T)` 表；加载时导出热力学状态表，运行时 C1 双三次 Hermite 插值 | 表中气体，且 `R` / 摩尔质量与工质匹配 |
| `helmholtz` | NIST He-4 Helmholtz | NIST He-4 Helmholtz 基本方程 | 仅 He |

输运模型独立于 EOS：`sutherland` 使用 Sutherland 黏度与固定 Pr；`nist-helium` 使用 He 专用 NIST 黏度与 Hands-Arp 导热，且仅适用于 He。

## 统一接口与单位

- `gas.state(T,rho)`：输入温度 `T` [K]、质量密度 `rho` [kg/m³]。
- `gas.density(T,p)` / `gas.stateTP(T,p)`：输入温度 [K]、压力 [Pa]。
- `gas.transport(T,p)` / `gas.transportAtDensity(T,rho)`：返回 `[mu,k,pr]`，单位为 Pa·s、W/(m·K)、无量纲。

所有 EOS 的稳定状态字段为 `Pressure`、`InternalEnergy`、`Enthalpy`、`Cv`、`Cp`、`PressureTemperatureDerivative`、`PressureDensityDerivative`、`SoundSpeed`、`Z`、`R`、`Density`，均为 SI、质量基准。`Entropy` 是 Tabular Gas 的可选扩展字段；不同 EOS 的绝对焓/熵参考态未统一，不能跨模型直接比较绝对值。

## Tabular Gas

Tabular Gas 接收已解析的 MAT 文件路径、SAGE 风格文本路径，或内存 struct。最小原始表字段如下：

```matlab
raw = struct( ...
    'TemperatureK',T, ...                  % nT-by-1，严格递增
    'SpecificVolumeM3PerKg',v, ...         % 1-by-nV，严格递增
    'Compressibility',Z, ...               % nT-by-nV
    'ReferenceCpJPerKgK',1039);            % 标量或 nT-by-1
gas = gasprop.create("N2","tabular",raw);
```

若表中没有 `InternalEnergyJPerKg` / `EntropyJPerKgK`，包会根据 SAGE 的 `Z(v,T)` 关系在**加载阶段**派生并检查它们；循环计算阶段只做短小的曲面求值和有界密度反演。表外推会被拒绝。本库不附带具体气体表数据。

## 验证（默认关闭）

验证代码单独放在 [validation](validation)，主包不会自动加入该目录、也不会自动运行验证。需要时明确执行：

```matlab
root = "D:\\VibeCodingWorkSpace\\GasPropertyResearch";
addpath(root)
addpath(fullfile(root,"validation"))
report = runGasPropertyValidation();
```

验证只检查纯物性契约：Ideal/RK/He 状态与反演、He 输运、Tabular 的 MAT/SAGE 文本加载和导数一致性；不包含 SIMPLE 工程的求解器或循环集成测试。
