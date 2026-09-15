classdef NISTHelium
    %HELIUMHELMHOLTZ NIST IR8474 normal-fluid He-4, independent of transport.
    properties (Constant)
        R=8.3144598/.004002602
    end
    methods
        function s=state(~,T,rho)
            gasprop.checkRange(T,[20,1500],'gasprop:HeliumEOSRange','Helium temperature [K]');
            gasprop.checkRange(rho,[realmin,500],'gasprop:HeliumEOSRange','Helium density [kg/m^3]');
            if gasprop.validation(),gaspropcheck.eosHelium('stateInput',T,rho);end
            s=gasprop.eos.internal.heliumState(T,rho);
            s.Density=rho+zeros(size(s.Pressure));
            gasprop.checkRange(s.Pressure,[realmin,2e9],'gasprop:HeliumEOSState','Helium pressure [Pa]');
            if gasprop.validation(),gaspropcheck.eosHelium('state',s);end
        end
        function rho=density(obj,T,p)
            gasprop.checkRange(T,[20,1500],'gasprop:HeliumEOSRange','Helium temperature [K]');
            gasprop.checkRange(p,[realmin,2e9],'gasprop:HeliumEOSPressure','Helium pressure [Pa]');
            if gasprop.validation(),gaspropcheck.eosHelium('densityInput',T,p);end
            rho=gasprop.eos.internal.invertDensity(obj,T,p,500);
        end
    end
end
