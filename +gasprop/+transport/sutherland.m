function [mu,k,pr]=sutherland(T,gas,cp)
%SUTHERLAND Dilute viscosity correlation with constant Prandtl number.
if gasprop.validation(),gaspropcheck.eosInputs(T);end
mu=gas.ReferenceViscosity.* ...
    (gas.ReferenceTemperature+gas.SutherlandTemperature)./ ...
    (T+gas.SutherlandTemperature).*(T/gas.ReferenceTemperature).^1.5;
if nargout>1
    if nargin<3, cp=gas.Cp; end
    mu=mu+zeros(size(T+cp));
    k=mu.*cp/gas.Prandtl;
    pr=gas.Prandtl+zeros(size(k));
end
end
