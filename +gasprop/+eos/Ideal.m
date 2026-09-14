classdef Ideal
    %IDEAL Calorically perfect gas, SI mass-specific properties.
    properties (SetAccess=private)
        R
        Cv
        Cp
    end
    methods
        function obj=Ideal(R,Cv,Cp)
            obj.R=R; obj.Cv=Cv; obj.Cp=Cp;
        end
        function s=state(obj,T,rho)
            gasprop.eos.validateInputs(T,rho);
            shape=T+rho; T=T+zeros(size(shape)); rho=rho+zeros(size(shape));
            p=rho.*obj.R.*T; cv=obj.Cv+zeros(size(T));
            s=struct('Pressure',p,'InternalEnergy',obj.Cv*T, ...
                'Enthalpy',obj.Cp*T,'Cv',cv,'Cp',obj.Cp+zeros(size(T)), ...
                'PressureTemperatureDerivative',rho*obj.R, ...
                'PressureDensityDerivative',obj.R*T, ...
                'SoundSpeed',sqrt(obj.Cp/obj.Cv*obj.R*T), ...
                'Z',ones(size(T)),'R',obj.R,'Density',rho);
        end
        function rho=density(obj,T,p)
            gasprop.eos.validateInputs(T,p);
            rho=p./(obj.R*T);
        end
    end
end
