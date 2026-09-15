function eosRK(action,obj,varargin)
%EOSRK Optional RK configuration, coverage and state diagnostics.
switch action
    case 'reference'
        data=varargin{1};
        if ~isstruct(data) || ~isscalar(data)
            error('gasprop:RKReferenceData','RK reference data must be a scalar structure.');
        end
        if isfield(data,'GasName') && ~strcmpi(string(data.GasName),obj.GasName)
            error('gasprop:RKReferenceData','Reference data gas must match %s.',obj.GasName);
        end
        fields={'GasConstant','CriticalTemperatureK','CriticalPressurePa', ...
            'ReferenceTemperatureK','MinimumTemperatureK','MaximumTemperatureK'};
        for k=1:numel(fields)
            if isfield(data,fields{k})
                validateattributes(data.(fields{k}),{'numeric'}, ...
                    {'real','finite','positive','scalar'});
            end
        end
        if isfield(data,'GasConstant') && abs(data.GasConstant-obj.R)>1e-10*obj.R
            error('gasprop:RKReferenceData', ...
                'Reference GasConstant must match the selected GasProperties.R.');
        end
        if isfield(data,'PressureRangePa')
            validateattributes(data.PressureRangePa,{'numeric'}, ...
                {'real','finite','positive','numel',2,'increasing'});
        end
        if ~all(isfield(data,{'TemperatureK','IdealCpJPerKgK'}))
            error('gasprop:RKReferenceData','Reference data requires TemperatureK and IdealCpJPerKgK.');
        end
        T=data.TemperatureK(:);cp=data.IdealCpJPerKgK(:);
        gaspropcheck.eosInputs(T,cp);
        if numel(T)<2 || numel(cp)~=numel(T) || any(diff(T)<=0) || any(cp<=obj.R)
            error('gasprop:RKReferenceData','Cp0 requires increasing T and matching Cp0>R values.');
        end
        referenceT=obj.ReferenceTemperatureK;
        if isfield(data,'ReferenceTemperatureK'),referenceT=data.ReferenceTemperatureK;end
        if referenceT<T(1) || referenceT>T(end)
            error('gasprop:RKReferenceData','Reference temperature lies outside Cp0 data.');
        end
    case 'domain'
        if obj.MinimumTemperatureK>=obj.MaximumTemperatureK || ...
                obj.PressureRangePa(1)>=obj.PressureRangePa(2)
            error('gasprop:RKReferenceData','Reference data leaves no usable gas-state domain.');
        end
    case 'stateInput'
        T=varargin{1};rho=varargin{2};
        gaspropcheck.eosInputs(T,rho);temperature(obj,T);
        if any(obj.B*rho(:)>=1)
            error('gasprop:RKRange','RK density requires 0<rho<1/b.');
        end
    case 'densityInput'
        T=varargin{1};p=varargin{2};
        gaspropcheck.eosInputs(T,p);temperature(obj,T);pressure(obj,p);
    case 'pressure'
        pressure(obj,varargin{1});
    case 'state'
        pr=varargin{1};cv=varargin{2};cp=varargin{3};c2=varargin{4};
        values=[pr(:);cv(:);cp(:);c2(:)];
        if ~isreal(values) || any(~isfinite(values) | values<=0)
            error('gasprop:RKState','Unstable or invalid RK state.');
        end
    case 'idealCp'
        T=varargin{1};range=varargin{2};
        gaspropcheck.eosInputs(T);temperature(obj,T);
        if ~isempty(range) && any(T(:)<range(1) | T(:)>range(2))
            error('gasprop:RKRange','Temperature lies outside the Cp0 reference table.');
        end
end
end

function temperature(obj,T)
if any(T(:)<obj.MinimumTemperatureK | T(:)>obj.MaximumTemperatureK)
    error('gasprop:RKRange','%s RK requires %.8g<=T<=%.8g K (gas/supercritical domain).', ...
        obj.GasName,obj.MinimumTemperatureK,obj.MaximumTemperatureK);
end
end

function pressure(obj,p)
tolerance=5e-10;
if any(~isfinite(p(:)) | p(:)<obj.PressureRangePa(1)*(1-tolerance) | ...
        p(:)>obj.PressureRangePa(2)*(1+tolerance))
    error('gasprop:RKRange','%s RK requires %.8g<=p<=%.8g Pa.', ...
        obj.GasName,obj.PressureRangePa(1),obj.PressureRangePa(2));
end
end
