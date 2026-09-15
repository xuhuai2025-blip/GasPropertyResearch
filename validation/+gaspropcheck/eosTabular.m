function eosTabular(action,table,varargin)
%EOSTABULAR Optional table-query coverage and state diagnostics.
switch action
    case 'stateInput'
        T=varargin{1};rho=varargin{2};
        gaspropcheck.eosInputs(T,rho);
        temperature(table,T);v=1./rho;
        if any(v(:)<table.SpecificVolumeM3PerKg(1) | v(:)>table.SpecificVolumeM3PerKg(end))
            error('gasprop:TabularRange','Tabular Gas state is outside the supplied v-T table range.');
        end
    case 'densityInput'
        T=varargin{1};p=varargin{2};
        gaspropcheck.eosInputs(T,p);temperature(table,T);pressure(table,p);
    case 'pressure'
        pressure(table,varargin{1});
    case 'bracket'
        p=varargin{1};lower=varargin{2};upper=varargin{3};
        if any(p<lower*(1-1e-12) | p>upper*(1+1e-12))
            error('gasprop:TabularPressureRange', ...
                'Requested T-p state is outside the pressure range bracketed by the table.');
        end
    case 'state'
        s=varargin{1};pressure(table,s.Pressure);
        values=[s.Pressure(:);s.Cv(:);s.Cp(:); ...
            s.PressureDensityDerivative(:);s.SoundSpeed(:)];
        if ~isreal(values) || any(~isfinite(values) | values<=0) || ...
                any(~isfinite(s.InternalEnergy(:)))
            error('gasprop:TabularState', ...
                'Tabular interpolation produced an invalid or unstable thermodynamic state.');
        end
end
end

function temperature(table,T)
if any(T(:)<table.TemperatureK(1) | T(:)>table.TemperatureK(end))
    error('gasprop:TabularRange','Tabular Gas state is outside the supplied v-T table range.');
end
end

function pressure(table,p)
if ~isempty(table.PressureRangePa)
    bounds=table.PressureRangePa;
    if any(~isfinite(p(:)) | p(:)<bounds(1)*(1-1e-9) | p(:)>bounds(2)*(1+1e-9))
        error('gasprop:TabularPressureRange', ...
            'Tabular Gas requires %.12g <= p <= %.12g Pa.',bounds(1),bounds(2));
    end
end
end
