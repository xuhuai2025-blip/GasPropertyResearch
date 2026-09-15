function tableBuild(stage,varargin)
%TABLEBUILD Optional raw-table checks and node stability checks.
switch stage
    case 'source'
        source(varargin{:});
    case 'state'
        state(varargin{:});
end
end

function source(raw,gas,sourceName)
if ~isstruct(raw) || ~isscalar(raw)
    error('gasprop:InvalidTabularData','Tabular Gas data must be one scalar struct.');
end
T=double(readField(raw,["TemperatureK","Temperature"],sourceName)); T=T(:);
v=double(readField(raw,["SpecificVolumeM3PerKg","SpecificVolume"],sourceName)); v=v(:).';
z=double(readField(raw,["Compressibility","Z"],sourceName));
matrix(z,numel(T),numel(v),'Compressibility',sourceName);
if numel(T)<4 || numel(v)<4 || ~isreal(T) || ~isreal(v) || ...
        any(~isfinite(T)) || any(~isfinite(v)) || any(T<=0) || any(v<=0) || ...
        any(diff(T)<=0) || any(diff(v)<=0) || any(z(:)<=0)
    error('gasprop:InvalidTabularData', ...
        'T, v and Z must be finite positive grids with at least 4 strictly ascending nodes.');
end
if isfield(raw,'GasConstant')
    checkGasConstant(raw.GasConstant,gas.R,sourceName);
elseif isfield(raw,'R')
    checkGasConstant(raw.R,gas.R,sourceName);
end
checkGasName(raw,gas,sourceName);
fields={'PressureTemperatureDerivativePaPerK','PressureDensityDerivativePaM3PerKg'};
if isfield(raw,'InternalEnergyJPerKg')
    fields=[fields,{'InternalEnergyJPerKg','CvJPerKgK','EntropyJPerKgK'}];
end
for i=1:numel(fields)
    if isfield(raw,fields{i})
        matrix(double(raw.(fields{i})),numel(T),numel(v),fields{i},sourceName);
    end
end
referenceCp(raw,gas,T,sourceName);
if isfield(raw,'PressureRangePa')
    bounds=double(raw.PressureRangePa(:).');
    validateattributes(bounds,{'numeric'},{'real','finite','positive','numel',2,'increasing'});
end
if isfield(raw,'TransportTable')
    tr=raw.TransportTable;
    t=double(tr.TemperatureK(:)); p=double(tr.PressurePa(:).');
    validateattributes(t,{'numeric'},{'real','finite','positive','increasing'});
    validateattributes(p,{'numeric'},{'real','finite','positive','increasing'});
    mu=double(tr.ViscosityPaS); k=double(tr.ThermalConductivityWPerMK);
    matrix(mu,numel(t),numel(p),'ViscosityPaS',sourceName);
    matrix(k,numel(t),numel(p),'ThermalConductivityWPerMK',sourceName);
    if any(mu(:)<=0 | k(:)<=0)
        error('gasprop:InvalidTabularData','Transport tables must be finite and positive.');
    end
end
end

function state(energy,entropy,energyT,p,T,v,z,zV,R)
if any(~isfinite(energy(:))) || any(~isfinite(entropy(:)))
    error('gasprop:InvalidTabularData', ...
        'Derived or supplied internal-energy/entropy tables must be finite.');
end
pRho=R.*T.*(z-v.*zV);
if any(p(:)<=0) || any(pRho(:)<=0) || any(energyT(:)<=0)
    error('gasprop:UnstableTabularData', ...
        ['Tabular Gas requires a mechanically and thermally stable single-phase ' ...
         'table: p>0, dp/drho>0, Cv>0 at every grid node.']);
end
end

function value=readField(raw,names,sourceName)
for name=names
    if isfield(raw,char(name)), value=raw.(char(name)); return, end
end
error('gasprop:InvalidTabularData', ...
    '%s is missing required field %s.',sourceName,strjoin(names,' or '));
end

function matrix(value,rowCount,columnCount,name,sourceName)
if ~isequal(size(value),[rowCount,columnCount]) || ~isreal(value) || any(~isfinite(value(:)))
    error('gasprop:InvalidTabularData', ...
        '%s in %s must be a real finite matrix matching its two grid axes.',name,sourceName);
end
end

function checkGasConstant(value,expected,sourceName)
value=double(value);
if ~(isscalar(value) && isfinite(value) && value>0) || abs(value-expected)>1e-10*expected
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
expected=strtrim(supplied);
if isscalar(expected) && strlength(expected)>0
    for entry=gasprop.gasCatalog()
        if any(lower(expected)==entry.Aliases)
            expected=entry.Name;
            break
        end
    end
end
if ~isscalar(expected) || (strlength(expected)>0 && ~strcmpi(expected,gas.Name))
    error('gasprop:TabularGasMismatch', ...
        '%s is for %s, but the operating case selected %s.',sourceName,expected,gas.Name);
end
if isfield(raw,'MolarMass')
    value=double(raw.MolarMass);
    if ~(isscalar(value) && isfinite(value) && value>0) || ...
            abs(value-gas.MolarMass)>1e-10*gas.MolarMass
        error('gasprop:TabularGasMismatch', ...
            '%s has a molar mass inconsistent with the selected Gas.',sourceName);
    end
end
end

function referenceCp(raw,gas,T,sourceName)
if isfield(raw,'ReferenceCpJPerKgK')
    cp0=double(raw.ReferenceCpJPerKgK);
elseif isfield(raw,'IdealCpJPerKgK')
    cp0=double(raw.IdealCpJPerKgK);
else
    cp0=gas.Cp;
    warning('gasprop:TabularCpFallback', ...
        ['%s does not supply cp0(T); Tabular Gas is using the selected ' ...
         'GasProperties constant Cp. Use a MAT table with ReferenceCpJPerKgK ' ...
         'for temperature-dependent cp0.'],sourceName);
end
if isscalar(cp0), cp0=cp0+zeros(size(T)); else, cp0=cp0(:); end
if numel(cp0)~=numel(T) || any(~isfinite(cp0)) || any(cp0<=gas.R)
    error('gasprop:InvalidTabularData', ...
        '%s ReferenceCpJPerKgK must be scalar or one positive value per T node, above R.',sourceName);
end
end
