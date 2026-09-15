classdef GasProperties
    %GASPROPERTIES 工质物性入口；EOS 与输运独立选择，He 默认 Helmholtz/NIST。
    % Cp/Cv/Gamma 为工质名录常量；RK 的 Cv0 在模型元数据中。
    % 实际状态比热必须通过 state/stateTP 获取。

    properties (SetAccess = private)
        Name
        MolarMass
        R
        Gamma
        Cv
        Cp
        ReferenceViscosity
        ReferenceTemperature
        SutherlandTemperature
        Prandtl
        EquationOfState = "ideal"
        TransportModel = "sutherland"
    end
    properties (Access=private)
        EOS
        ReferenceData = []
    end

    methods
        function obj = GasProperties(name, molarMass, gasConstant, gamma, ...
                referenceViscosity, referenceTemperature, ...
                sutherlandTemperature, prandtl)
            if gasprop.validation()
                gaspropcheck.gasConstants(name,[molarMass,gasConstant,gamma, ...
                    referenceViscosity,referenceTemperature,sutherlandTemperature,prandtl]);
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
            obj.EOS = gasprop.eos.IdealGas(obj.R,obj.Cv,obj.Cp);
        end

        function value = viscosity(obj, temperature, pressure)
            %VISCOSITY [Pa*s]. NIST and tabular models require pressure.
            if any(obj.TransportModel==["nist-helium","tabular"])
                if nargin<3
                    error('gasprop:PressureRequired','The selected viscosity model requires temperature and pressure.');
                end
                rho=obj.density(temperature,pressure);
                if obj.TransportModel=="nist-helium"
                    [value,~]=gasprop.transport.heliumTransport(temperature,rho);
                else
                    [value,~]=obj.EOS.transport(temperature,rho,pressure);
                end
                return
            end
            if obj.TransportModel=="reference"
                [value,~]=gasprop.transport.reference(temperature,obj.ReferenceData);
            else
                value=gasprop.transport.sutherland(temperature,obj);
            end
        end

        function obj=withTransportModel(obj,name)
            name=gasprop.GasProperties.normalizeTransportModel(name);
            if name=="nist-helium" && obj.Name~="Helium"
                error('gasprop:TransportGasCombination','NIST helium transport supports He only.');
            end
            if name=="tabular" && (obj.EquationOfState~="tabular" || ~obj.EOS.StateTable.HasTransport)
                error('gasprop:MissingTransportData','Tabular transport requires a selected table with transport data.');
            end
            if name=="reference"
                if isempty(obj.ReferenceData)
                    obj.ReferenceData=gasprop.table.rkReference(obj.Name);
                end
                if isempty(obj.ReferenceData)
                    error('gasprop:MissingTransportData','No bundled reference transport data for %s.',obj.Name);
                end
            end
            obj.TransportModel=name;
        end

        function [mu,k,pr]=transport(obj,T,p)
            %TRANSPORT SI transport at T,p; reference Pr uses Cp0, others use state Cp.
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
            % Optional transport envelope check, controlled by gasprop.validation.
            if gasprop.validation() && obj.TransportModel=="nist-helium"
                gaspropcheck.transport('helium',T,rho);
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
                    model=gasprop.eos.IdealGas(obj.R,obj.Cv,obj.Cp);
                case "rk"
                    if isempty(tabularData), tabularData=gasprop.table.rkReference(obj.Name); end
                    model=gasprop.eos.RKGas(obj,tabularData);
                    obj.ReferenceData=tabularData;
                case "tabular"
                    if isempty(tabularData)
                        catalog=gasprop.gasCatalog();
                        entry=find(string({catalog.Name})==obj.Name,1);
                        if ~isempty(entry)
                            root=fileparts(fileparts(mfilename('fullpath')));
                            tabularData=fullfile(root,'data',catalog(entry).TabularFile);
                        else
                            error('gasprop:MissingTabularData','Tabular Gas requires an explicit table for %s.',obj.Name);
                        end
                    end
                    model=gasprop.eos.TabularGas(obj,tabularData);
                case "helmholtz"
                    if obj.Name~="Helium"
                        error('gasprop:EOSGasCombination', ...
                            'NIST He-4 Helmholtz supports He only.');
                    end
                    model=gasprop.eos.NISTHelium();
            end
            obj.EOS=model; obj.EquationOfState=name;
            if obj.TransportModel=="tabular" && name~="tabular"
                obj.TransportModel="sutherland";
            end
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
            if any(obj.EquationOfState==["tabular","rk"])
                value=obj.EOS.metadata();
            end
        end

        function value=idealCp(obj,T)
            %IDEALCP Ideal reference Cp; distinct from state-dependent Cp.
            if obj.EquationOfState=="rk"
                value=obj.EOS.idealCp(T);
            else
                if gasprop.validation(), gaspropcheck.eosInputs(T); end
                value=obj.Cp+zeros(size(T));
            end
        end

        function obj=withDefaultTransport(obj)
            %WITHDEFAULTTRANSPORT Select transport accompanying bundled data.
            if obj.EquationOfState=="tabular" && obj.EOS.StateTable.HasTransport
                obj=obj.withTransportModel("tabular");
            elseif obj.EquationOfState=="rk" && ~isempty(obj.ReferenceData)
                obj=obj.withTransportModel("reference");
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
            if ~isscalar(name) || ismissing(name) || ~any(name==["sutherland","nist-helium","tabular","reference"])
                error('gasprop:TransportModel','Transport must be sutherland, nist-helium, tabular or reference.');
            end
        end

        function label=transportModelDescription(name)
            name=gasprop.GasProperties.normalizeTransportModel(name);
            if name=="sutherland"
                label="Sutherland 黏度 + 固定 Pr";
            elseif name=="nist-helium"
                label="NIST 氦气黏度 + Hands-Arp 导热";
            elseif name=="tabular"
                label="REFPROP 表格黏度与导热系数";
            else
                label="零密度参考物性随温度插值";
            end
        end

        function gas = fromName(name)
            inputName=lower(strtrim(string(name)));
            canonical="";
            if isscalar(inputName) && ~ismissing(inputName)
                for entry=gasprop.gasCatalog()
                    if any(inputName==entry.Aliases), canonical=entry.Name; break, end
                end
            end
            switch canonical
                case "Helium"
                    gas = gasprop.GasProperties( ...
                        "Helium", 0.004002602, gasprop.eos.NISTHelium.R, 5/3, ...
                        18.85e-6, 273, 80, 0.71);
                    gas=gas.withEquationOfState('helmholtz').withTransportModel('nist-helium');
                case "Nitrogen"
                    gas = gasprop.GasProperties( ...
                        "Nitrogen", 28.0134e-3, 296.8, 1.4, ...
                        16.63e-6, 273, 107, 0.71);
                case "Argon"
                    gas = gasprop.GasProperties( ...
                        "Argon", 39.948e-3, 208.1, 1.67, ...
                        21.25e-6, 273, 148, 0.67);
                case "Hydrogen"
                    gas = gasprop.GasProperties( ...
                        "Hydrogen", 2.01588e-3, 8.314472/2.01588e-3, 1.4, ...
                        8.35e-6, 273, 84.4, 0.71);
                case "Air"
                    gas = gasprop.GasProperties( ...
                        "Air", 28.97e-3, 287.0, 1.4, ...
                        17.08e-6, 273, 112.0, 0.71);
                case "CarbonDioxide"
                    % REFPROP 9.1 CO2.FLD: ideal Cp and dilute mu/k at 300 K.
                    % Sutherland S matches dilute viscosity at 300 and 305 K;
                    % it is a local ideal-model approximation, not a full-range fit.
                    R=8.31451/.0440098;
                    cp0=845.84601038795631;
                    mu0=1.5013923982513168e-5;
                    k0=.016747259697437032;
                    ratio=(1.5255769724942452e-5/mu0)/(305/300)^1.5;
                    sutherland=(300-305*ratio)/(ratio-1);
                    gas=gasprop.GasProperties("CarbonDioxide",.0440098,R, ...
                        cp0/(cp0-R),mu0,300,sutherland,mu0*cp0/k0);
                otherwise
                    error('gasprop:UnsupportedGas', ...
                        'Unsupported gas name. Use He, H2, N2, Ar, Air or CO2.');
            end
        end
    end

    methods (Access=private)
        function [mu,k,pr]=transportForState(obj,T,s)
            if obj.TransportModel=="nist-helium"
                [mu,k]=gasprop.transport.heliumTransport(T,s.Density);
                pr=mu.*s.Cp./k;
            elseif obj.TransportModel=="tabular"
                [mu,k]=obj.EOS.transport(T,s.Density,s.Pressure);
                pr=mu.*s.Cp./k;
            elseif obj.TransportModel=="reference"
                [mu,k,cp0]=gasprop.transport.reference(T,obj.ReferenceData);
                mu=mu+zeros(size(s.Pressure)); k=k+zeros(size(s.Pressure));
                pr=mu.*cp0./k;
            else
                [mu,k,pr]=gasprop.transport.sutherland(T,obj,s.Cp);
            end
        end
    end
end
