classdef GasProperties
    %GASPROPERTIES 工质物性入口；EOS 与输运独立选择，He 默认 Helmholtz/NIST。
    % Cp/Cv/Gamma 为理想参考常量，真实气体物性必须通过 state/stateTP 获取。

    properties (SetAccess = private)
        Name (1,1) string
        MolarMass (1,1) double
        R (1,1) double
        Gamma (1,1) double
        Cv (1,1) double
        Cp (1,1) double
        ReferenceViscosity (1,1) double
        ReferenceTemperature (1,1) double
        SutherlandTemperature (1,1) double
        Prandtl (1,1) double
        EquationOfState (1,1) string = "ideal"
        TransportModel (1,1) string = "sutherland"
    end
    properties (Access=private)
        EOS
    end

    methods
        function obj = GasProperties(name, molarMass, gasConstant, gamma, ...
                referenceViscosity, referenceTemperature, ...
                sutherlandTemperature, prandtl)
            values = [molarMass, gasConstant, gamma, referenceViscosity, ...
                referenceTemperature, sutherlandTemperature, prandtl];
            validateattributes(values, {'numeric'}, ...
                {'real', 'finite', 'positive', 'vector'});
            if gamma <= 1
                error('gasprop:InvalidGamma', ...
                    'Specific heat ratio must be greater than one.');
            end
            obj.Name = string(name);
            obj.MolarMass = molarMass;
            obj.R = gasConstant;
            obj.Gamma = gamma;
            obj.Cv = gasConstant / (gamma - 1);
            obj.Cp = gamma * obj.Cv;
            obj.ReferenceViscosity = referenceViscosity;
            obj.ReferenceTemperature = referenceTemperature;
            obj.SutherlandTemperature = sutherlandTemperature;
            obj.Prandtl = prandtl;
            obj.EOS = gasprop.eos.Ideal(obj.R,obj.Cv,obj.Cp);
        end

        function value = viscosity(obj, temperature, pressure)
            %VISCOSITY [Pa*s]. NIST requires pressure; no silent dilute fallback.
            if obj.TransportModel=="nist-helium"
                if nargin<3
                    error('gasprop:PressureRequired','NIST helium viscosity requires T,p; use heliumTransport(T,0) explicitly for the dilute limit.');
                end
                [value,~]=gasprop.heliumTransport(temperature,obj.density(temperature,pressure));
                return
            end
            validateattributes(temperature, {'numeric'}, ...
                {'real', 'finite', 'positive'});
            value = obj.ReferenceViscosity .* ...
                (obj.ReferenceTemperature + obj.SutherlandTemperature) ./ ...
                (temperature + obj.SutherlandTemperature) .* ...
                (temperature ./ obj.ReferenceTemperature).^1.5;
        end

        function obj=withTransportModel(obj,name)
            name=gasprop.GasProperties.normalizeTransportModel(name);
            if name=="nist-helium" && obj.Name~="Helium"
                error('gasprop:TransportGasCombination','NIST helium transport supports He only.');
            end
            obj.TransportModel=name;
        end

        function [mu,k,pr]=transport(obj,T,p)
            %TRANSPORT SI transport at T,p using density and Cp from selected EOS.
            if nargin<3
                if obj.TransportModel~="sutherland" || obj.EquationOfState~="ideal"
                    error('gasprop:PressureRequired','State-dependent transport requires pressure.');
                end
                mu=obj.viscosity(T);k=mu*obj.Cp/obj.Prandtl;
                pr=obj.Prandtl+zeros(size(mu));return
            end
            s=obj.stateTP(T,p);
            [mu,k,pr]=obj.transportForState(T,s);
        end

        function [mu,k,pr]=transportAtDensity(obj,T,rho)
            %TRANSPORTATDENSITY Use known local mass density without TP inversion.
            s=obj.state(T,rho);
            [mu,k,pr]=obj.transportForState(T,s);
        end

        function validateTransportState(obj,T,rho)
            % Checks the combined operating envelope even at integrator stages.
            if obj.TransportModel=="nist-helium"
                gasprop.validateHeliumTransportInputs(T,rho);
            end
        end

        function value = specificEnthalpy(obj, temperature, pressure)
            if obj.EquationOfState=="ideal"
                value = obj.Cp .* temperature;
            else
                if nargin<3
                    error('gasprop:PressureRequired', ...
                        'Real-gas enthalpy requires temperature and pressure.');
                end
                s=obj.stateTP(temperature,pressure); value=s.Enthalpy;
            end
        end

        function obj=withEquationOfState(obj,name,tabularData)
            %WITHEQUATIONOFSTATE Select one of the four gas-model families.
            %   Tabular Gas receives a SAGE-style Z(v,T) state-table source.
            if nargin<3, tabularData=[]; end
            name=gasprop.GasProperties.normalizeEquationOfState(name);
            switch name
                case "ideal"
                    model=gasprop.eos.Ideal(obj.R,obj.Cv,obj.Cp);
                case "rk"
                    model=gasprop.eos.RedlichKwong(obj);
                case "tabular"
                    if isempty(tabularData)
                        error('gasprop:MissingTabularData', ...
                            'Tabular Gas requires a TableDataFile containing Z(v,T).');
                    end
                    model=gasprop.eos.TabularGas(obj,tabularData);
                case "helmholtz"
                    if obj.Name~="Helium"
                        error('gasprop:EOSGasCombination', ...
                            'NIST He-4 Helmholtz supports He only.');
                    end
                    model=gasprop.eos.HeliumHelmholtz();
            end
            obj.EOS=model; obj.EquationOfState=name;
        end

        function s=state(obj,T,rho)
            s=obj.EOS.state(T,rho);
        end

        function rho=density(obj,T,p)
            rho=obj.EOS.density(T,p);
        end

        function s=stateTP(obj,T,p)
            s=obj.state(T,obj.density(T,p));
        end

        function value=equationOfStateMetadata(obj)
            value=struct('Name',obj.EquationOfState, ...
                'Description',gasprop.GasProperties.equationOfStateDescription( ...
                    obj.EquationOfState));
            if obj.EquationOfState=="tabular"
                value=obj.EOS.metadata();
            end
        end
    end

    methods (Static)
        function name=normalizeEquationOfState(name)
            name=lower(strtrim(string(name)));
            if ~isscalar(name) || ismissing(name)
                error('gasprop:EquationOfState', ...
                    'EOS must be ideal, rk, tabular or helmholtz.');
            end
            switch name
                case "ideal gas"
                    name="ideal";
                case "rk gas"
                    name="rk";
                case "tabular gas"
                    name="tabular";
                case {"nist he-4 helmholtz","nist-he4-helmholtz"}
                    name="helmholtz";
            end
            if ~any(name==["ideal","rk","tabular","helmholtz"])
                error('gasprop:EquationOfState', ...
                    'EOS must be ideal, rk, tabular or helmholtz.');
            end
        end

        function label=equationOfStateDescription(name)
            name=gasprop.GasProperties.normalizeEquationOfState(name);
            switch name
                case "ideal"
                    label="Ideal Gas (p = rho R T)";
                case "rk"
                    label="RK Gas (Redlich-Kwong)";
                case "tabular"
                    label="Tabular Gas (SAGE Z(v,T) bicubic)";
                case "helmholtz"
                    label="NIST He-4 Helmholtz (He only)";
            end
        end

        function name=normalizeTransportModel(name)
            name=lower(strtrim(string(name)));
            if ~isscalar(name) || ismissing(name) || ~any(name==["sutherland","nist-helium"])
                error('gasprop:TransportModel','Transport must be sutherland or nist-helium.');
            end
        end

        function label=transportModelDescription(name)
            name=gasprop.GasProperties.normalizeTransportModel(name);
            if name=="sutherland"
                label="Sutherland 黏度 + 固定 Pr";
            else
                label="NIST 氦气黏度 + Hands-Arp 导热";
            end
        end

        function gas = fromName(name)
            switch upper(string(name))
                case "HE"
                    gas = gasprop.GasProperties( ...
                        "Helium", 0.004002602, gasprop.eos.HeliumHelmholtz.R, 5/3, ...
                        18.85e-6, 273, 80, 0.71);
                    gas=gas.withEquationOfState('helmholtz').withTransportModel('nist-helium');
                case "N2"
                    gas = gasprop.GasProperties( ...
                        "Nitrogen", 28.0134e-3, 296.8, 1.4, ...
                        16.63e-6, 273, 107, 0.71);
                case "AR"
                    gas = gasprop.GasProperties( ...
                        "Argon", 39.948e-3, 208.1, 1.67, ...
                        21.25e-6, 273, 148, 0.67);
                case "H2"
                    gas = gasprop.GasProperties( ...
                        "Hydrogen", 2.016e-3, 4157.2, 1.4, ...
                        8.35e-6, 273, 84.4, 0.71);
                case "AIR"
                    gas = gasprop.GasProperties( ...
                        "Air", 28.97e-3, 287.0, 1.4, ...
                        17.08e-6, 273, 112.0, 0.71);
                otherwise
                    error('gasprop:UnsupportedGas', ...
                        ['Unsupported gas "%s". Use He, H2, Air, N2, ' ...
                         'or Ar.'], name);
            end
        end
    end

    methods (Access=private)
        function [mu,k,pr]=transportForState(obj,T,s)
            if obj.TransportModel=="nist-helium"
                [mu,k]=gasprop.heliumTransport(T,s.Density);
                pr=mu.*s.Cp./k;
            else
                mu=obj.viscosity(T)+zeros(size(s.Pressure));
                k=mu.*s.Cp/obj.Prandtl;
                pr=obj.Prandtl+zeros(size(k));
            end
        end
    end
end
