function table=build(raw,gas,sourceName)
%BUILD Compile a raw T-v table once; optional checks live in validation/.
if nargin<3, sourceName="in-memory table"; end
if gasprop.validation(), gaspropcheck.tableBuild('source',raw,gas,sourceName); end
T=double(readField(raw,["TemperatureK","Temperature"])); T=T(:);
v=double(readField(raw,["SpecificVolumeM3PerKg","SpecificVolume"])); v=v(:).';
z=double(readField(raw,["Compressibility","Z"]));
zT=gasprop.table.internal.differentiate(z,T,1);
zV=gasprop.table.internal.differentiate(z,v,2);
if isfield(raw,'PressureTemperatureDerivativePaPerK')
    pT=double(raw.PressureTemperatureDerivativePaPerK);
    zT=(pT.*v/gas.R-z)./T;
end
if isfield(raw,'PressureDensityDerivativePaM3PerKg')
    pRho=double(raw.PressureDensityDerivativePaM3PerKg);
    zV=(z-pRho./(gas.R*T))./v;
end
zSurface=gasprop.table.internal.buildSurface(T,v,z,zT,zV);
[cp0,cp0Source]=referenceCp(raw,gas,T);
hasSuppliedEnergy=isfield(raw,'InternalEnergyJPerKg');
hasSuppliedEntropy=isfield(raw,'EntropyJPerKgK');
if hasSuppliedEnergy
    energy=double(raw.InternalEnergyJPerKg);
    if isfield(raw,'CvJPerKgK')
        energyT=double(raw.CvJPerKgK);
    else
        energyT=gasprop.table.internal.differentiate(energy,T,1);
    end
    energyV=gas.R.*T.^2.*zT./v;
    if hasSuppliedEntropy
        entropy=double(raw.EntropyJPerKgK);
    else
        entropy=gasprop.table.internal.deriveEntropy( ...
            T,v,z,gas.R,energyT(:,end),zSurface, ...
            gasprop.table.internal.buildSurface( ...
                T,v,energy,energyT,energyV));
    end
else
    [energy,entropy,energyV]=gasprop.table.internal.deriveStateTables( ...
        T,v,z,gas.R,cp0,zT,zSurface);
    energyT=gasprop.table.internal.differentiate(energy,T,1);
end
rho=1./v;
p=gas.R.*T.*z./v;
if gasprop.validation()
    gaspropcheck.tableBuild('state',energy,entropy,energyT,p,T,v,z,zV,gas.R);
end
energySurface=gasprop.table.internal.buildSurface( ...
    T,v,energy,energyT,energyV);
entropySurface=gasprop.table.internal.buildSurface( ...
    T,v,entropy,energyT./T,(energyV+p)./T);
table=struct('SchemaVersion',"2.0",'Kind',"SAGE-TabularGas", ...
    'Source',sourceMetadata(raw,sourceName),'SourceFile',string(sourceName), ...
    'GasName',gas.Name,'GasConstant',gas.R, ...
    'TemperatureK',T,'SpecificVolumeM3PerKg',v, ...
    'Compressibility',z,'CompressibilitySurface',zSurface, ...
    'InternalEnergyJPerKg',energy,'InternalEnergySurface',energySurface, ...
    'EntropyJPerKgK',entropy,'ReferenceCpJPerKgK',cp0, ...
    'ReferenceCpSource',cp0Source, ...
    'EntropySurface',entropySurface, ...
    'MinimumDensityKgPerM3',min(rho),'MaximumDensityKgPerM3',max(rho), ...
    'MinimumPressurePa',min(p(:)),'MaximumPressurePa',max(p(:)));
table.PressureRangePa=[];
if isfield(raw,'PressureRangePa')
    table.PressureRangePa=double(raw.PressureRangePa(:).');
end
if isfield(raw,'SourceMetadata'), table.SourceMetadata=raw.SourceMetadata; end
table.HasTransport=isfield(raw,'TransportTable');
if table.HasTransport
    transport=raw.TransportTable;
    transport.TemperatureK=double(transport.TemperatureK(:));
    transport.PressurePa=double(transport.PressurePa(:).');
    transport.ViscosityPaS=double(transport.ViscosityPaS);
    transport.ThermalConductivityWPerMK=double(transport.ThermalConductivityWPerMK);
    table.TransportTable=transport;
end
end

function value=readField(raw,names)
for name=names
    if isfield(raw,char(name)), value=raw.(char(name)); return, end
end
value=raw.(char(names(1)));
end

function value=sourceMetadata(raw,sourceName)
if isfield(raw,'Source') && strlength(strtrim(string(raw.Source)))>0
    value=string(raw.Source);
else
    value=string(sourceName);
end
end

function [cp0,source]=referenceCp(raw,gas,T)
if isfield(raw,'ReferenceCpJPerKgK')
    cp0=double(raw.ReferenceCpJPerKgK);
    source="table ReferenceCpJPerKgK";
elseif isfield(raw,'IdealCpJPerKgK')
    cp0=double(raw.IdealCpJPerKgK);
    source="table IdealCpJPerKgK";
else
    cp0=gas.Cp;
    source="selected-gas constant Cp (calorically-perfect fallback)";
end
if isscalar(cp0), cp0=cp0+zeros(size(T)); else, cp0=cp0(:); end
end