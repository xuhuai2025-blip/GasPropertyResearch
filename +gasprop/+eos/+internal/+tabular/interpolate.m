function [value,temperatureDerivative,volumeDerivative]=interpolate(surface,v,T)
%INTERPOLATE Evaluate one C1 bicubic surface; optional checks are at the facade.
if nargout<2
    value=gasprop.eos.internal.tabular.evaluateSurface(surface,T,v);
elseif nargout<3
    [value,temperatureDerivative]=gasprop.eos.internal.tabular.evaluateSurface(surface,T,v);
else
    [value,temperatureDerivative,volumeDerivative]= ...
        gasprop.eos.internal.tabular.evaluateSurface(surface,T,v);
end
end
