function transport(action,varargin)
%TRANSPORT Optional transport input, coverage and state diagnostics.
switch action
    case 'helium'
        T=varargin{1};rho=varargin{2};
        validateattributes(T,{'numeric'},{'real','nonempty'});
        validateattributes(rho,{'numeric'},{'real','nonempty'});
        if any(~isfinite(T(:)) | T(:)<20 | T(:)>830) || ...
                any(~isfinite(rho(:)) | rho(:)<0 | rho(:)>160)
            error('gasprop:HeliumTransportRange', ...
                'Helium transport requires 20 <= T <= 830 K and 0 <= rho <= 160 kg/m^3.');
        end
    case 'reference'
        T=varargin{1};grid=varargin{2};
        gaspropcheck.eosInputs(T);
        if any(T(:)<grid(1) | T(:)>grid(end))
            error('gasprop:ReferenceTemperatureRange', ...
                'Reference properties require %.12g <= T <= %.12g K.',grid(1),grid(end));
        end
    case 'tabular'
        table=varargin{1};T=varargin{2};p=varargin{3};
        gaspropcheck.eosInputs(T,p);
        gaspropcheck.eosTabular('pressure',table,p);
        bounds=table.TransportTable.PressurePa([1,end]);
        temperatures=table.TransportTable.TemperatureK([1,end]);
        if any(p(:)<bounds(1)*(1-1e-9) | p(:)>bounds(2)*(1+1e-9)) || ...
                any(T(:)<temperatures(1) | T(:)>temperatures(2))
            error('gasprop:TabularTransportRange', ...
                'Transport table requires %.12g..%.12g K and %.12g..%.12g Pa.', ...
                temperatures(1),temperatures(2),bounds(1),bounds(2));
        end
    case 'state'
        mu=varargin{1};k=varargin{2};values=[mu(:);k(:)];
        if ~isreal(values) || any(~isfinite(values) | values<=0)
            error('gasprop:TabularTransportState','Interpolated transport must be finite and positive.');
        end
end
end
