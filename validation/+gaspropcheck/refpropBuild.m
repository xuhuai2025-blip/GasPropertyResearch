function refpropBuild(kind,varargin)
%REFPROPBUILD Offline sampling checks, enabled by buildRefpropTables.
switch kind
    case 'options'
        options=varargin{1};
        for name=["TemperatureCount","VolumeCount","PressureCount","ReferenceCount"]
            validateattributes(options.(name),{'numeric'}, ...
                {'scalar','real','finite','integer','>=',4});
        end
    case 'melting'
        [fluid,tmelt,ierr,herr]=varargin{:};
        if ierr~=0 || ~isfinite(tmelt) || tmelt<=0
            error('refpropbuild:MissingPhaseBoundary', ...
                'No usable melting boundary for %s at 20 MPa: %s', ...
                fluid.Name,strtrim(char(herr(:)')));
        end
    case 'meltingRange'
        [tmelt,bounds]=varargin{:};
        if numel(bounds)~=2 || tmelt<bounds(1) || tmelt>bounds(2)
            error('refpropbuild:MeltingRange', ...
                '20 MPa melting point is outside the declared model range.');
        end
    case 'domain'
        [minimum,maximum,fluid]=varargin{:};
        if minimum>=maximum
            error('refpropbuild:UnsupportedRange','No supported gas range for %s.',fluid.Name);
        end
    case 'padding'
        [info,T,rho,p,tm]=varargin{:};
        if ~isfinite(tm) || tm<info.MeltingModelTemperatureRangeK(1) || ...
                tm>info.MeltingModelTemperatureRangeK(2) || T<=tm || ...
                p>info.Limits.EOS.MaximumPressurePa || ...
                rho/(info.REFPROPMolarMassKgPerMol*1000)>info.Limits.EOS.MaximumMolarDensityMolPerL
            error('refpropbuild:UnsafePadding', ...
                'T-v padding is outside supported fluid/EOS bounds at T=%g K.',T);
        end
    case 'node'
        [T,rho,p,u,s,cv,pT,drdp,q]=varargin{:};
        if (q>=0 && q<=1) || ~all(isfinite([p,u,s,cv,pT,drdp])) || any([p,cv,drdp]<=0)
            error('refpropbuild:InvalidNode','Unsupported state at T=%g K, rho=%g kg/m3.',T,rho);
        end
    case 'transportPressure'
        [info,p]=varargin{:};
        if p>min(info.Limits.ETA.MaximumPressurePa,info.Limits.TCX.MaximumPressurePa)
            error('refpropbuild:TransportRange','Requested pressure exceeds a transport model limit.');
        end
    case 'transportNode'
        [info,T,p,mu,k,rho,q]=varargin{:};
        if (q>=0 && q<=1) || ~all(isfinite([mu,k])) || any([mu,k]<=0) || ...
                rho/(info.REFPROPMolarMassKgPerMol*1000)>min(info.Limits.ETA.MaximumMolarDensityMolPerL,info.Limits.TCX.MaximumMolarDensityMolPerL)
            error('refpropbuild:TransportRange','Unsupported transport node at T=%g K p=%g Pa.',T,p);
        end
    case 'reference'
        values=vertcat(varargin{:});
        if any(~isfinite(values) | values<=0)
            error('refpropbuild:InvalidReference','REFPROP returned invalid ideal reference properties.');
        end
end
end
