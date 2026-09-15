function [mu,k,cp0]=reference(T,data)
%REFERENCE Interpolate offline zero-density REFPROP transport and ideal Cp.
grid=data.TemperatureK;
if gasprop.validation(),gaspropcheck.transport('reference',T,grid);end
mu=interp1(grid,data.ViscosityPaS,T,'pchip');
k=interp1(grid,data.ThermalConductivityWPerMK,T,'pchip');
cp0=interp1(grid,data.IdealCpJPerKgK,T,'pchip');
end
