classdef NISTHelium
    %HELIUMHELMHOLTZ NIST IR8474 normal-fluid He-4, independent of transport.
    properties (Constant)
        R=8.3144598/.004002602
    end
    methods
        function s=state(~,T,rho)
            if gasprop.validation(),gaspropcheck.eosHelium('stateInput',T,rho);end
            s=gasprop.eos.internal.heliumState(T,rho);
            s.Density=rho+zeros(size(s.Pressure));
            if gasprop.validation(),gaspropcheck.eosHelium('state',s);end
        end
        function rho=density(obj,T,p)
            if gasprop.validation(),gaspropcheck.eosHelium('densityInput',T,p);end
            rho=gasprop.eos.internal.invertDensity(obj,T,p,500);
        end
    end
end
