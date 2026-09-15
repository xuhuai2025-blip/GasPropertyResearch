classdef RKGas
    %RKGAS Redlich-Kwong EOS with constant reference ideal Cv.
    % See docs/sage-rk-comparison.md. Only stable gas/supercritical states.
    % A Cp0(T) reference table sets Cv at T0; idealCp exposes the separate
    % dilute-property curve without changing the energy derivatives.
    properties (SetAccess=private)
        GasName
        R
        IdealCv
        CriticalTemperature
        CriticalPressure
        EffectiveCriticalTemperature
        A
        B
        C1
        C2
        ReferenceTemperatureK = 300
        MinimumTemperatureK
        MaximumTemperatureK = 800
        PressureRangePa = [1e5 20e6]
        ParameterSource
        IsPseudoPure
    end
    properties (Access=private)
        CpCurve = []
        CpTemperatureRange = []
    end
    methods
        function obj=RKGas(gas,referenceData)
            if nargin<2, referenceData=[]; end
            parameters=gasprop.eos.internal.rkParameters(gas.Name);
            obj.GasName=string(gas.Name); obj.R=gas.R;
            obj.IdealCv=gas.Cv;
            obj.C1=parameters.C1; obj.C2=parameters.C2;
            obj.CriticalTemperature=parameters.CriticalTemperatureK;
            obj.CriticalPressure=parameters.CriticalPressurePa;
            obj.ParameterSource=parameters.Source;
            obj.IsPseudoPure=parameters.IsPseudoPure;
            dataMinimum=10;
            if ~isempty(referenceData)
                if gasprop.validation(),gaspropcheck.eosRK('reference',obj,referenceData);end
                if isfield(referenceData,'CriticalTemperatureK')
                    obj.CriticalTemperature=referenceData.CriticalTemperatureK;
                end
                if isfield(referenceData,'CriticalPressurePa')
                    obj.CriticalPressure=referenceData.CriticalPressurePa;
                end
                if isfield(referenceData,'ReferenceTemperatureK')
                    obj.ReferenceTemperatureK=referenceData.ReferenceTemperatureK;
                end
                if isfield(referenceData,'MinimumTemperatureK')
                    dataMinimum=max(dataMinimum,referenceData.MinimumTemperatureK);
                end
                if isfield(referenceData,'MaximumTemperatureK')
                    obj.MaximumTemperatureK=min(800,referenceData.MaximumTemperatureK);
                end
                if isfield(referenceData,'PressureRangePa')
                    range=referenceData.PressureRangePa;
                    obj.PressureRangePa=[max(1e5,range(1)), min(20e6,range(2))];
                end
                if isfield(referenceData,'TemperatureK') && isfield(referenceData,'IdealCpJPerKgK')
                    T=referenceData.TemperatureK(:); cp=referenceData.IdealCpJPerKgK(:);
                    obj.CpCurve=pchip(T,cp);
                    obj.CpTemperatureRange=[T(1),T(end)];
                    obj.IdealCv=ppval(obj.CpCurve,obj.ReferenceTemperatureK)-obj.R;
                    dataMinimum=max(dataMinimum,T(1));
                    obj.MaximumTemperatureK=min(obj.MaximumTemperatureK,T(end));
                end
            end
            obj.A=obj.C1*obj.R^2*obj.CriticalTemperature^2.5/obj.CriticalPressure;
            obj.B=obj.C2*obj.R*obj.CriticalTemperature/obj.CriticalPressure;
            % Fitting a,b moves the cubic EOS critical point. In particular,
            % the Sage H2 fit has a higher critical T than the physical fluid.
            obj.EffectiveCriticalTemperature=(obj.A/(obj.R*obj.B)* ...
                .08664034996495772/.4274802335403414)^(2/3);
            obj.MinimumTemperatureK=max(dataMinimum, ...
                1.001*max(obj.CriticalTemperature,obj.EffectiveCriticalTemperature));
            if gasprop.validation(),gaspropcheck.eosRK('domain',obj);end
        end

        function s=state(obj,T,rho)
            gasprop.checkRange(T,[obj.MinimumTemperatureK,obj.MaximumTemperatureK],'gasprop:RKRange','RK temperature [K]');
            if gasprop.validation(),gaspropcheck.eosRK('stateInput',obj,T,rho);end
            shape=T+rho; T=T+zeros(size(shape)); rho=rho+zeros(size(shape));
            [p,pT,pr]=obj.pressureTerms(T,rho);
            gasprop.checkRange(p,obj.PressureRangePa,'gasprop:RKRange','RK pressure [Pa]',5e-10);
            if gasprop.validation(),gaspropcheck.eosRK('pressure',obj,p);end
            a=obj.A; b=obj.B; br=b*rho;
            u=obj.IdealCv*T-1.5*a/b./sqrt(T).*log1p(br);
            cv=obj.IdealCv+.75*a/b./T.^1.5.*log1p(br);
            cp=cv+T.*pT.^2./(rho.^2.*pr);
            c2=pr+T.*pT.^2./(rho.^2.*cv);
            if gasprop.validation(),gaspropcheck.eosRK('state',obj,pr,cv,cp,c2);end
            entropy=obj.IdealCv*log(T)+obj.R*log(1./rho-b) ...
                -.5*a/b./T.^1.5.*log1p(br);
            s=struct('Pressure',p,'InternalEnergy',u,'Enthalpy',u+p./rho, ...
                'Entropy',entropy,'Cv',cv,'Cp',cp,'PressureTemperatureDerivative',pT, ...
                'PressureDensityDerivative',pr,'SoundSpeed',sqrt(c2), ...
                'Z',p./(rho*obj.R.*T),'R',obj.R,'Density',rho);
        end

        function rho=density(obj,T,p)
            % Algebraic gas root of Sage equation (18.18), then roundoff polish.
            gasprop.checkRange(T,[obj.MinimumTemperatureK,obj.MaximumTemperatureK],'gasprop:RKRange','RK temperature [K]');
            gasprop.checkRange(p,obj.PressureRangePa,'gasprop:RKRange','RK pressure [Pa]',5e-10);
            if gasprop.validation(),gaspropcheck.eosRK('densityInput',obj,T,p);end
            shape=T+p; T=T+zeros(size(shape)); p=p+zeros(size(shape));
            ar=obj.A*p./(obj.R^2*T.^2.5); br=obj.B*p./(obj.R*T);
            coefficient=ar-br.^2-br;
            depressedP=coefficient-1/3;
            depressedQ=-2/27+coefficient/3-ar.*br;
            discriminant=(depressedQ/2).^2+(depressedP/3).^3;
            z=zeros(size(p)); one=discriminant>=0;
            root=sqrt(max(discriminant(one),0));
            u=-depressedQ(one)/2+root; v=-depressedQ(one)/2-root;
            z(one)=sign(u).*abs(u).^(1/3)+sign(v).*abs(v).^(1/3)+1/3;
            three=~one;
            if any(three(:))
                radius=sqrt(-depressedP(three)/3);
                argument=-depressedQ(three)./(2*radius.^3);
                z(three)=2*radius.*cos(acos(max(-1,min(1,argument)))/3)+1/3;
            end
            rho=p./(z*obj.R.*T);
            for k=1:2
                [calculated,~,derivative]=obj.pressureTerms(T,rho);
                rho=rho-(calculated-p)./derivative;
            end
            [calculated,~,derivative]=obj.pressureTerms(T,rho);
            if any(~isfinite(rho(:)) | derivative(:)<=0 | abs(calculated(:)-p(:))./p(:)>2e-10)
                error('gasprop:EOSInversion','RK density inversion did not converge to a stable gas root.');
            end
        end

        function [p,pT,pr]=pressureTerms(obj,T,rho)
            % Internal EOS primitives; request-domain checks belong to callers.
            br=obj.B*rho;
            attraction=obj.A*rho.^2./(sqrt(T).*(1+br));
            p=obj.R*T.*rho./(1-br)-attraction;
            pT=obj.R*rho./(1-br)+attraction./(2*T);
            pr=obj.R*T./(1-br).^2-obj.A./sqrt(T).*rho.*(2+br)./(1+br).^2;
        end

        function value=idealCp(obj,T)
            % Temperature-dependent dilute Cp for property/transport queries.
            gasprop.checkRange(T,[obj.MinimumTemperatureK,obj.MaximumTemperatureK],'gasprop:RKRange','RK reference temperature [K]');
            if gasprop.validation(),gaspropcheck.eosRK('idealCp',obj,T,obj.CpTemperatureRange);end
            if isempty(obj.CpCurve)
                value=obj.IdealCv+obj.R+zeros(size(T));
            else
                value=ppval(obj.CpCurve,T);
            end
        end

        function value=metadata(obj)
            value=struct('Name',"rk",'Description',"RK Gas (Redlich-Kwong)", ...
                'GasName',obj.GasName,'C1',obj.C1,'C2',obj.C2, ...
                'ParameterSource',obj.ParameterSource,'IsPseudoPure',obj.IsPseudoPure, ...
                'CriticalTemperatureK',obj.CriticalTemperature, ...
                'CriticalPressurePa',obj.CriticalPressure, ...
                'EffectiveCriticalTemperatureK',obj.EffectiveCriticalTemperature, ...
                'TemperatureRangeK',[obj.MinimumTemperatureK,obj.MaximumTemperatureK], ...
                'PressureRangePa',obj.PressureRangePa, ...
                'ReferenceTemperatureK',obj.ReferenceTemperatureK, ...
                'IdealCvJPerKgK',obj.IdealCv,'HasReferenceCpTable',~isempty(obj.CpCurve));
        end
    end
end
