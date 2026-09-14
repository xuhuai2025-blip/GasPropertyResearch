function table=normalizeStateTable(raw,gas,sourceName)
%NORMALIZESTATETABLE Validate one immutable SAGE-style tabular-gas table.
if ~isstruct(raw) || ~isscalar(raw)
    error('gasprop:InvalidTabularData', ...
        'Tabular Gas data must be one scalar struct.');
end
T=column(readField(raw,["TemperatureK","Temperature"],sourceName));
v=row(readField(raw,["SpecificVolumeM3PerKg","SpecificVolume"],sourceName));
z=matrix(readField(raw,["Compressibility","Z"],sourceName), ...
    numel(T),numel(v),'Compressibility',sourceName);
if numel(T)<4 || numel(v)<4 || any(~isfinite(T)) || any(~isfinite(v)) || ...
        any(T<=0) || any(v<=0) || any(diff(T)<=0) || any(diff(v)<=0) || ...
        any(~isfinite(z(:))) || any(z(:)<=0)
    error('gasprop:InvalidTabularData', ...
        'T, v and Z must be finite positive grids with at least 4 strictly ascending nodes.');
end
if isfield(raw,'GasConstant')
    checkGasConstant(raw.GasConstant,gas.R,sourceName);
elseif isfield(raw,'R')
    checkGasConstant(raw.R,gas.R,sourceName);
end
checkGasName(raw,gas,sourceName);
zT=gasprop.eos.tabular.differentiate(z,T,1);
zV=gasprop.eos.tabular.differentiate(z,v,2);
zSurface=gasprop.eos.tabular.buildSurface(T,v,z,zT,zV);
[cp0,cp0Source]=referenceCp(raw,gas,T,sourceName);
hasSuppliedEnergy=isfield(raw,'InternalEnergyJPerKg');
hasSuppliedEntropy=isfield(raw,'EntropyJPerKgK');
if hasSuppliedEnergy
    energy=matrix(raw.InternalEnergyJPerKg,numel(T),numel(v), ...
        'InternalEnergyJPerKg',sourceName);
    energyT=gasprop.eos.tabular.differentiate(energy,T,1);
    energyV=gasprop.eos.tabular.differentiate(energy,v,2);
    if hasSuppliedEntropy
        entropy=matrix(raw.EntropyJPerKgK,numel(T),numel(v), ...
            'EntropyJPerKgK',sourceName);
    else
        entropy=gasprop.eos.tabular.deriveEntropy( ...
            T,v,z,gas.R,energyT(:,end),zSurface, ...
            gasprop.eos.tabular.buildSurface( ...
                T,v,energy,energyT,energyV));
    end
else
    [energy,entropy,energyV]=gasprop.eos.tabular.deriveStateTables( ...
        T,v,z,gas.R,cp0,zT,zSurface);
    energyT=gasprop.eos.tabular.differentiate(energy,T,1);
end
if any(~isfinite(energy(:))) || any(~isfinite(entropy(:)))
    error('gasprop:InvalidTabularData', ...
        'Derived or supplied internal-energy/entropy tables must be finite.');
end
rho=1./v;
p=gas.R.*T.*z./v;
pRho=gas.R.*T.*(z-v.*zV);
if any(p(:)<=0) || any(pRho(:)<=0) || any(energyT(:)<=0)
    error('gasprop:UnstableTabularData', ...
        ['Tabular Gas requires a mechanically and thermally stable single-phase ' ...
         'table: p>0, dp/drho>0, Cv>0 at every grid node.']);
end
energySurface=gasprop.eos.tabular.buildSurface( ...
    T,v,energy,energyT,energyV);
entropySurface=gasprop.eos.tabular.buildSurface( ...
    T,v,entropy,energyT./T,(energyV+p)./T);
gasprop.eos.tabular.validateSurfaceState( ...
    zSurface,energySurface,gas.R);
if hasSuppliedEnergy
    validateSuppliedThermodynamics(zSurface,energySurface,entropySurface, ...
        gas.R,hasSuppliedEntropy,sourceName);
end
table=struct('SchemaVersion',"1.0",'Kind',"SAGE-TabularGas", ...
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
end

function value=readField(raw,names,sourceName)
for name=names
    if isfield(raw,char(name)), value=raw.(char(name)); return, end
end
error('gasprop:InvalidTabularData', ...
    '%s is missing required field %s.',sourceName,strjoin(names,' or '));
end

function value=column(value)
value=double(value(:));
end

function value=row(value)
value=double(value(:).');
end

function value=matrix(value,rowCount,columnCount,name,sourceName)
value=double(value);
if ~isequal(size(value),[rowCount,columnCount])
    error('gasprop:InvalidTabularData', ...
        '%s in %s must be numel(TemperatureK)-by-numel(SpecificVolumeM3PerKg).', ...
        name,sourceName);
end
end

function checkGasConstant(value,expected,sourceName)
value=double(value);
if ~(isscalar(value) && isfinite(value) && value>0) || ...
        abs(value-expected)>1e-10*expected
    error('gasprop:TabularGasMismatch', ...
        '%s has a gas constant inconsistent with the selected Gas.',sourceName);
end
end

function checkGasName(raw,gas,sourceName)
if isfield(raw,'GasName')
    supplied=string(raw.GasName);
elseif isfield(raw,'Gas')
    supplied=string(raw.Gas);
else
    supplied="";
end
if strlength(strtrim(supplied))==0
    checkMolarMass(raw,gas,sourceName);
    return
end
switch lower(strtrim(supplied))
    case {"he","he4","he-4","helium"}, expected="Helium";
    case {"h2","hydrogen"}, expected="Hydrogen";
    case {"n2","nitrogen"}, expected="Nitrogen";
    case {"ar","argon"}, expected="Argon";
    case "air", expected="Air";
    otherwise
        error('gasprop:TabularGasMismatch', ...
            '%s names an unsupported table gas "%s".',sourceName,supplied);
end
if expected~=gas.Name
    error('gasprop:TabularGasMismatch', ...
        '%s is for %s, but the operating case selected %s.', ...
        sourceName,expected,gas.Name);
end
checkMolarMass(raw,gas,sourceName);
end

function checkMolarMass(raw,gas,sourceName)
if isfield(raw,'MolarMass')
    value=double(raw.MolarMass);
    if ~(isscalar(value) && isfinite(value) && value>0) || ...
            abs(value-gas.MolarMass)>1e-10*gas.MolarMass
        error('gasprop:TabularGasMismatch', ...
            '%s has a molar mass inconsistent with the selected Gas.',sourceName);
    end
end
end

function value=sourceMetadata(raw,sourceName)
if isfield(raw,'Source') && strlength(strtrim(string(raw.Source)))>0
    value=string(raw.Source);
else
    value=string(sourceName);
end
end

function validateSuppliedThermodynamics(zSurface,energySurface,entropySurface, ...
        R,hasSuppliedEntropy,sourceName)
fractions=[0.211324865405187,0.5,0.788675134594813];
tLower=zSurface.TemperatureK(1:end-1);
vLower=zSurface.SpecificVolumeM3PerKg(1:end-1).';
T=sort(reshape(tLower+(zSurface.TemperatureK(2:end)-tLower).* ...
    fractions,[],1));
v=sort(reshape(vLower+(zSurface.SpecificVolumeM3PerKg(2:end).'-vLower).* ...
    fractions,[],1)).';
[v,T]=meshgrid(v,T);
[z,zT]=gasprop.eos.tabular.interpolate(zSurface,v,T);
[~,energyT,energyV]=gasprop.eos.tabular.interpolate( ...
    energySurface,v,T);
pressure=R.*T.*z./v;
if relativeMismatch(energyV,R.*T.^2.*zT./v)>0.05
    error('gasprop:InconsistentTabularEnergy', ...
        ['%s supplied epsilon(v,T) is inconsistent with Z(v,T): ' ...
         'epsilon_v must equal R*T^2*Z_T/v.'],sourceName);
end
if hasSuppliedEntropy
    [~,entropyT,entropyV]=gasprop.eos.tabular.interpolate( ...
        entropySurface,v,T);
    if relativeMismatch(entropyT,energyT./T)>0.05 || ...
            relativeMismatch(entropyV,(energyV+pressure)./T)>0.05
        error('gasprop:InconsistentTabularEntropy', ...
            ['%s supplied s(v,T) is inconsistent with epsilon(v,T) and ' ...
             'Z(v,T).'],sourceName);
    end
end
end

function value=relativeMismatch(actual,expected)
scale=max(max(1,abs(actual)),abs(expected));
value=max(abs(actual-expected)./scale,[],'all');
end

function [cp0,source]=referenceCp(raw,gas,T,sourceName)
if isfield(raw,'ReferenceCpJPerKgK')
    cp0=double(raw.ReferenceCpJPerKgK);
    source="table ReferenceCpJPerKgK";
elseif isfield(raw,'IdealCpJPerKgK')
    cp0=double(raw.IdealCpJPerKgK);
    source="table IdealCpJPerKgK";
else
    cp0=gas.Cp;
    source="selected-gas constant Cp (calorically-perfect fallback)";
    warning('gasprop:TabularCpFallback', ...
        ['%s does not supply cp0(T); Tabular Gas is using the selected ' ...
         'GasProperties constant Cp. Use a MAT table with ReferenceCpJPerKgK ' ...
         'for temperature-dependent cp0.'],sourceName);
end
if isscalar(cp0), cp0=cp0+zeros(size(T)); else, cp0=cp0(:); end
if numel(cp0)~=numel(T) || any(~isfinite(cp0)) || any(cp0<=gas.R)
    error('gasprop:InvalidTabularData', ...
        '%s ReferenceCpJPerKgK must be scalar or one positive value per T node, above R.', ...
        sourceName);
end
end
