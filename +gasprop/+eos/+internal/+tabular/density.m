function rho=density(table,R,T,p)
%DENSITY SAGE-style pressure bracket followed by safeguarded Newton steps.
shape=size(T+p); %#ok<ELARLOG> Deliberate scalar expansion to common state shape.
T=reshape(T+zeros(shape),1,[]); p=reshape(p+zeros(shape),1,[]);
rhoLower=ones(size(T))./table.SpecificVolumeM3PerKg(end);
rhoUpper=ones(size(T))./table.SpecificVolumeM3PerKg(1);
if gasprop.validation()
    lowerPressure=gasprop.eos.internal.tabular.pressure(table,R,T,rhoLower);
    upperPressure=gasprop.eos.internal.tabular.pressure(table,R,T,rhoUpper);
    gaspropcheck.eosTabular('bracket',table,p,lowerPressure,upperPressure);
end
rho=min(max(p./(R.*T),rhoLower),rhoUpper);
lo=rhoLower; hi=rhoUpper;
for iteration=1:60
    [currentPressure,~,pressureDerivative]=gasprop.eos.internal.tabular.pressure(table,R,T,rho);
    residual=currentPressure-p;
    if max(abs(residual)./p)<1e-10
        rho=reshape(rho,shape); return
    end
    hi(residual>0)=rho(residual>0);
    lo(residual<=0)=rho(residual<=0);
    candidate=rho-residual./pressureDerivative;
    outside=candidate<=lo | candidate>=hi | ~isfinite(candidate);
    candidate(outside)=(lo(outside)+hi(outside))/2;
    done=abs(residual)./p<1e-10;
    candidate(done)=rho(done);
    rho=candidate;
end
error('gasprop:TabularDensityInversion', ...
    'Tabular Gas T-p to density inversion did not converge.');
end
