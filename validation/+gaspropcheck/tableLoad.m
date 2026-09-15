function tableLoad(stage,value,gas)
%TABLELOAD Optional source and compiled-table format checks.
switch stage
    case 'source'
        if ~(isstruct(value) || isstring(value) || ischar(value))
            error('gasprop:InvalidTabularData', ...
                'Tabular Gas source must be a table struct or a .mat/.txt file.');
        end
    case 'text'
        if size(value,1)<5 || size(value,2)<5
            error('gasprop:InvalidTabularData', ...
                'SAGE Z(v,T) text table needs one header row/column and at least 4-by-4 data.');
        end
    case 'compiled'
        compiled(value,gas);
end
end

function compiled(table,gas)
required={'Kind','GasName','GasConstant','TemperatureK','SpecificVolumeM3PerKg', ...
    'Source','SourceFile','ReferenceCpSource', ...
    'CompressibilitySurface','InternalEnergySurface','EntropySurface', ...
    'MinimumDensityKgPerM3','MaximumDensityKgPerM3'};
if ~any(string(table.SchemaVersion)==["1.0","2.0"]) || ...
        ~all(isfield(table,required)) || string(table.Kind)~="SAGE-TabularGas"
    error('gasprop:InvalidCompiledTable','Unsupported or incomplete compiled gas table.');
end
if string(table.GasName)~=gas.Name || ~isscalar(table.GasConstant) || ...
        ~isfinite(table.GasConstant) || abs(table.GasConstant-gas.R)>1e-10*gas.R
    error('gasprop:TabularGasMismatch','Compiled table does not match the selected gas.');
end
T=table.TemperatureK; v=table.SpecificVolumeM3PerKg;
validateattributes(T,{'numeric'},{'column','real','finite','positive','increasing'});
validateattributes(v,{'numeric'},{'row','real','finite','positive','increasing'});
if numel(T)<4 || numel(v)<4
    error('gasprop:InvalidCompiledTable','Compiled table requires at least four nodes per axis.');
end
fields={'CompressibilitySurface','InternalEnergySurface','EntropySurface'};
parts={'Value','TemperatureDerivative','VolumeDerivative','CrossDerivative'};
for i=1:numel(fields)
    surface=table.(fields{i});
    if ~isstruct(surface) || ~all(isfield(surface,[parts,{'TemperatureK','SpecificVolumeM3PerKg'}])) || ...
            ~isequal(surface.TemperatureK,T) || ~isequal(surface.SpecificVolumeM3PerKg,v)
        error('gasprop:InvalidCompiledTable','Compiled surface axes do not match the table.');
    end
    for j=1:numel(parts)
        values=surface.(parts{j});
        if ~isequal(size(values),[numel(T),numel(v)]) || ~isreal(values) || any(~isfinite(values(:)))
            error('gasprop:InvalidCompiledTable','Invalid shape or value in compiled surface.');
        end
    end
end
if ~isempty(table.PressureRangePa)
    validateattributes(table.PressureRangePa,{'numeric'}, ...
        {'real','finite','positive','numel',2,'increasing'});
end
if table.HasTransport
    if ~isfield(table,'TransportTable')
        error('gasprop:InvalidCompiledTable','Missing T-p transport table.');
    end
    tr=table.TransportTable;
    validateattributes(tr.TemperatureK,{'numeric'},{'column','real','finite','positive','increasing'});
    validateattributes(tr.PressurePa,{'numeric'},{'row','real','finite','positive','increasing'});
    expected=[numel(tr.TemperatureK),numel(tr.PressurePa)];
    if ~isequal(size(tr.ViscosityPaS),expected) || ~isequal(size(tr.ThermalConductivityWPerMK),expected) || ...
            any(~isfinite(tr.ViscosityPaS(:)) | tr.ViscosityPaS(:)<=0 | ...
            ~isfinite(tr.ThermalConductivityWPerMK(:)) | tr.ThermalConductivityWPerMK(:)<=0)
        error('gasprop:InvalidCompiledTable','Invalid T-p transport values.');
    end
end
end
