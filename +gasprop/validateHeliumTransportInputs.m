function validateHeliumTransportInputs(T,rho)
%VALIDATEHELIUMTRANSPORTINPUTS Independent transport envelope, SI units.
validateattributes(T,{'numeric'},{'real','nonempty'});
validateattributes(rho,{'numeric'},{'real','nonempty'});
if any(~isfinite(T(:)) | T(:)<20 | T(:)>830) || ...
        any(~isfinite(rho(:)) | rho(:)<0 | rho(:)>160)
    error('gasprop:HeliumTransportRange', ...
        'Helium transport requires 20 <= T <= 830 K and 0 <= rho <= 160 kg/m^3.');
end
end
