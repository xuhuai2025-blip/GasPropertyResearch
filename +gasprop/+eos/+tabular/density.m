function rho=density(table,R,T,p)
%DENSITY SAGE-style pressure bracket followed by safeguarded Newton steps.
shape=size(T+p); %#ok<ELARLOG> Deliberate scalar expansion to common state shape.
T=reshape(T+zeros(shape),1,[]); p=reshape(p+zeros(shape),1,[]);
rhoLower=ones(size(T))./table.SpecificVolumeM3PerKg(end);
rhoUpper=ones(size(T))./table.SpecificVolumeM3PerKg(1);
lowerState=gasprop.eos.tabular.state(table,R,T,rhoLower);
upperState=gasprop.eos.tabular.state(table,R,T,rhoUpper);
if any(p<lowerState.Pressure | p>upperState.Pressure)
    error('gasprop:TabularPressureRange', ...
        'Requested T-p state is outside the pressure range bracketed by the table.');
end
rho=min(max(p./(R.*T),rhoLower),rhoUpper);
lo=rhoLower; hi=rhoUpper;
for iteration=1:60
    current=gasprop.eos.tabular.state(table,R,T,rho);
    residual=current.Pressure-p;
    if max(abs(residual)./p)<1e-10
        rho=reshape(rho,shape); return
    end
    hi(residual>0)=rho(residual>0);
    lo(residual<=0)=rho(residual<=0);
    candidate=rho-residual./current.PressureDensityDerivative;
    outside=candidate<=lo | candidate>=hi | ~isfinite(candidate);
    candidate(outside)=(lo(outside)+hi(outside))/2;
    rho=candidate;
end
error('gasprop:TabularDensityInversion', ...
    'Tabular Gas T-p to density inversion did not converge.');
end
