classdef HeliumHelmholtz
    %HELIUMHELMHOLTZ NIST IR8474 normal-fluid He-4, independent of transport.
    properties (Constant)
        R=8.3144598/.004002602
    end
    methods
        function s=state(~,T,rho)
            gasprop.eos.validateInputs(T,rho);
            s=gasprop.eos.heliumState(T,rho);
            s.Density=rho+zeros(size(s.Pressure));
        end
        function rho=density(obj,T,p)
            if any(p(:)>2e9)
                error('gasprop:HeliumEOSPressure','He pressure must not exceed 2 GPa.');
            end
            rho=gasprop.eos.invertDensity(obj,T,p,500);
        end
    end
end
