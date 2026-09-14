classdef RedlichKwong
    %REDLICHKWONG Standard RK with thermodynamically consistent departure terms.
    % p=R*T/(v-b)-a/(sqrt(T)*v*(v+b)); v in m3/kg.
    % u^res=-3*a/(2*b*sqrt(T))*log(1+b*rho).
    % Reference ideal Cv is inherited from the existing gas definition.
    % Single phase only: T>Tc, no liquid/metastable root selection.
    properties (SetAccess=private)
        R
        IdealCv
        CriticalTemperature
        CriticalPressure
        A
        B
    end
    methods
        function obj=RedlichKwong(gas)
            % Rounded critical points from CoolProp fluid pages (2026-09-09):
            % https://coolprop.org/fluid_properties/fluids/{name}.html
            switch gas.Name
                case "Helium", Tc=5.195300013635951; pc=228322.7892147868;
                case "Hydrogen", Tc=33.1443326883113; pc=1296357.6060553084;
                case "Nitrogen", Tc=126.192; pc=3395800;
                case "Argon", Tc=150.687; pc=4863000;
                otherwise
                    error('gasprop:EOSGasCombination', ...
                        'RK supports pure He, H2, N2 and Ar; Air has no configured mixing model.');
            end
            obj.R=gas.R; obj.IdealCv=gas.Cv;
            obj.CriticalTemperature=Tc; obj.CriticalPressure=pc;
            obj.A=.4274802335403414*gas.R^2*Tc^2.5/pc;
            obj.B=.08664034996495772*gas.R*Tc/pc;
        end
        function s=state(obj,T,rho)
            gasprop.eos.validateInputs(T,rho);
            shape=T+rho; T=T+zeros(size(shape)); rho=rho+zeros(size(shape));
            if any(T(:)<=obj.CriticalTemperature | T(:)<20 | T(:)>1500) || any(obj.B*rho(:)>=1)
                error('gasprop:RKRange', ...
                    'RK debug range: 20<=T<=1500 K, T>Tc and 0<rho<1/b. No two-phase model.');
            end
            a=obj.A; b=obj.B; R=obj.R; br=b*rho;
            attraction=a*rho.^2./(sqrt(T).*(1+br));
            p=R*T.*rho./(1-br)-attraction;
            pT=R*rho./(1-br)+attraction./(2*T);
            pr=R*T./(1-br).^2-a./sqrt(T).*rho.*(2+br)./(1+br).^2;
            u=obj.IdealCv*T-1.5*a/b./sqrt(T).*log1p(br);
            cv=obj.IdealCv+.75*a/b./T.^1.5.*log1p(br);
            cp=cv+T.*pT.^2./(rho.^2.*pr);
            c2=pr+T.*pT.^2./(rho.^2.*cv);
            if any(~isfinite(p(:)) | p(:)<=0 | pr(:)<=0 | cv(:)<=0 | ~isfinite(cp(:)))
                error('gasprop:RKState','Unstable or invalid RK state.');
            end
            s=struct('Pressure',p,'InternalEnergy',u,'Enthalpy',u+p./rho, ...
                'Cv',cv,'Cp',cp,'PressureTemperatureDerivative',pT, ...
                'PressureDensityDerivative',pr,'SoundSpeed',sqrt(c2), ...
                'Z',p./(rho*R.*T),'R',R,'Density',rho);
        end
        function rho=density(obj,T,p)
            rho=gasprop.eos.invertDensity(obj,T,p,(1-1e-10)/obj.B);
        end
    end
end
