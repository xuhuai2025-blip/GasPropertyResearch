function rho=invertDensity(eos,T,p,maxDensity)
%INVERTDENSITY Safeguarded Newton iteration for a stable single-phase EOS.
gasprop.eos.validateInputs(T,p);
shape=T+p; T=T+zeros(size(shape)); p=p+zeros(size(shape));
rho=min(p./(eos.R*T),.8*maxDensity);
lo=zeros(size(rho)); hi=maxDensity+zeros(size(rho));
for k=1:80
    try
        s=eos.state(T,rho);
    catch exception
        % At high pressure an ideal-density trial can exceed the He EOS
        % pressure limit even when the requested state is supported. Reduce
        % the trial, never clip the requested physical state or its output.
        if ~strcmp(exception.identifier,'gasprop:HeliumEOSState')
            rethrow(exception)
        end
        if ~isscalar(rho)
            rho=arrayfun(@(temperature,pressure)gasprop.eos.invertDensity( ...
                eos,temperature,pressure,maxDensity),T,p);
            return
        end
        hi=rho; rho=(lo+hi)/2; continue
    end
    err=s.Pressure-p;
    if max(abs(err./p),[],'all')<1e-11, return; end
    if any((maxDensity-rho(:))/maxDensity<1e-12 & err(:)<0)
        error('gasprop:EOSDensityRange', ...
            'Requested T,p requires density above the supported EOS range.');
    end
    hi(err>0)=rho(err>0); lo(err<=0)=rho(err<=0);
    step=err./s.PressureDensityDerivative;
    next=rho-max(min(step,.5*rho),-.5*rho);
    outside=next<=lo | next>=hi | ~isfinite(next);
    next(outside)=.5*(lo(outside)+hi(outside));
    done=abs(err./p)<1e-11;next(done)=rho(done);
    rho=next;
end
error('gasprop:EOSInversion','EOS density inversion did not converge.');
end
