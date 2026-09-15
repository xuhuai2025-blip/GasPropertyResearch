function rho=invertDensity(eos,T,p,maxDensity)
%INVERTDENSITY Safeguarded Newton iteration for the NIST helium EOS.
shape=T+p; T=T+zeros(size(shape)); p=p+zeros(size(shape));
rho=min(p./(eos.R*T),.8*maxDensity);
lo=zeros(size(rho)); hi=maxDensity+zeros(size(rho));
for k=1:80
    s=gasprop.eos.internal.heliumState(T,rho);
    % Shrink unsuitable trial states directly, independently of diagnostics.
    unsuitable=~isfinite(s.Pressure) | s.Pressure<=0 | s.Pressure>2e9 | ...
        ~isfinite(s.PressureDensityDerivative) | s.PressureDensityDerivative<=0;
    if any(unsuitable(:))
        hi(unsuitable)=rho(unsuitable);
        rho(unsuitable)=(lo(unsuitable)+hi(unsuitable))/2;
        continue
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
