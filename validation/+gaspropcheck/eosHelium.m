function eosHelium(action,varargin)
%EOSHELIUM Optional NIST He EOS input, coverage and state diagnostics.
switch action
    case 'stateInput'
        T=varargin{1};rho=varargin{2};
        gaspropcheck.eosInputs(T,rho);
        if any(T(:)<20 | T(:)>1500) || any(rho(:)>500)
            error('gasprop:HeliumEOSRange', ...
                'He EOS wrapper requires 20..1500 K and 0 < rho <= 500 kg/m3.');
        end
    case 'densityInput'
        T=varargin{1};p=varargin{2};
        gaspropcheck.eosInputs(T,p);
        if any(T(:)<20 | T(:)>1500)
            error('gasprop:HeliumEOSRange','He EOS wrapper requires 20..1500 K.');
        end
        if any(p(:)>2e9)
            error('gasprop:HeliumEOSPressure','He pressure must not exceed 2 GPa.');
        end
    case 'state'
        s=varargin{1};
        values=[s.Pressure(:);s.Cv(:);s.Cp(:); ...
            s.PressureDensityDerivative(:);s.SoundSpeed(:)];
        if ~isreal(values) || any(~isfinite(values) | values<=0) || any(s.Pressure(:)>2e9)
            error('gasprop:HeliumEOSState', ...
                'Invalid or out-of-pressure-range normal-fluid helium state.');
        end
end
end
