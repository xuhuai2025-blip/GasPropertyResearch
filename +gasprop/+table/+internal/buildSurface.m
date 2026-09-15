function surface=buildSurface(T,v,values,temperatureDerivative,volumeDerivative)
%BUILDSURFACE Cache one C1 bicubic-Hermite surface and its node gradients.
%   Gradients come from cubic splines at load time.  Evaluation always
%   returns a value and derivatives from this same local polynomial.
if nargin<4 || isempty(temperatureDerivative)
    temperatureDerivative=gasprop.table.internal.differentiate(values,T,1);
end
if nargin<5 || isempty(volumeDerivative)
    volumeDerivative=gasprop.table.internal.differentiate(values,v,2);
end
if gasprop.validation()
    gaspropcheck.tableSurface('surface',values,temperatureDerivative,volumeDerivative);
end
crossDerivative=(gasprop.table.internal.differentiate( ...
    temperatureDerivative,v,2)+gasprop.table.internal.differentiate( ...
    volumeDerivative,T,1))/2;
surface=struct('TemperatureK',T,'SpecificVolumeM3PerKg',v, ...
    'Value',values,'TemperatureDerivative',temperatureDerivative, ...
    'VolumeDerivative',volumeDerivative,'CrossDerivative',crossDerivative);
end
