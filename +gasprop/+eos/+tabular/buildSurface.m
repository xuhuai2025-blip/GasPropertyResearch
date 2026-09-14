function surface=buildSurface(T,v,values,temperatureDerivative,volumeDerivative)
%BUILDSURFACE Cache one C1 bicubic-Hermite surface and its node gradients.
%   Gradients come from cubic splines at load time.  Evaluation always
%   returns a value and derivatives from this same local polynomial.
if nargin<4 || isempty(temperatureDerivative)
    temperatureDerivative=gasprop.eos.tabular.differentiate(values,T,1);
end
if nargin<5 || isempty(volumeDerivative)
    volumeDerivative=gasprop.eos.tabular.differentiate(values,v,2);
end
if ~isequal(size(values),size(temperatureDerivative),size(volumeDerivative))
    error('gasprop:InvalidTabularData', ...
        'A tabular surface and its node gradients must have the same shape.');
end
crossDerivative=(gasprop.eos.tabular.differentiate( ...
    temperatureDerivative,v,2)+gasprop.eos.tabular.differentiate( ...
    volumeDerivative,T,1))/2;
surface=struct('TemperatureK',T,'SpecificVolumeM3PerKg',v, ...
    'Value',values,'TemperatureDerivative',temperatureDerivative, ...
    'VolumeDerivative',volumeDerivative,'CrossDerivative',crossDerivative);
end
