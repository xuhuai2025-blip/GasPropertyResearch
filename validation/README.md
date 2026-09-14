# Optional validation

此目录不在 `+gasprop` 包内，也不会由库自动加入 MATLAB path。仅在需要手动回归时执行：

```matlab
root = "D:\\VibeCodingWorkSpace\\GasPropertyResearch";
addpath(root)
addpath(fullfile(root,"validation"))
report = runGasPropertyValidation();
```

`runGasPropertyValidation` 只在系统临时目录创建小型解析表，不写入核心目录，也不调用 SIMPLE 求解器。`buildTabularGasTable` 是可选建表工具：它把 SAGE 风格 `Z(v,T)` 原始表编译成普通 struct；该 struct 可被 `gasprop.create(...,"tabular",table)` 或 `withEquationOfState` 直接使用。
