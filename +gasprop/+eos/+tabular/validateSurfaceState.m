function validateSurfaceState(zSurface,energySurface,R)
%VALIDATESURFACESTATE Reject spline-cell instability before solver runtime.
fractions=[0.211324865405187,0.5,0.788675134594813];
tLower=zSurface.TemperatureK(1:end-1);
vLower=zSurface.SpecificVolumeM3PerKg(1:end-1).';
T=sort(reshape(tLower+(zSurface.TemperatureK(2:end)-tLower).* ...
    fractions,[],1));
v=sort(reshape(vLower+(zSurface.SpecificVolumeM3PerKg(2:end).'-vLower).* ...
    fractions,[],1)).';
[v,T]=meshgrid(v,T);
[z,zT,zV]=gasprop.eos.tabular.interpolate(zSurface,v,T);
[~,cv]=gasprop.eos.tabular.interpolate(energySurface,v,T);
rho=1./v;
p=R.*T.*z./v;
pT=R.*(z+T.*zT)./v;
pRho=R.*T.*(z-v.*zV);
cp=cv+T.*pT.^2./(rho.^2.*pRho);
soundSquared=pRho+T.*pT.^2./(rho.^2.*cv);
if any(~isfinite(p),'all') || any(~isfinite(cv),'all') || ...
        any(~isfinite(pRho),'all') || any(~isfinite(cp),'all') || ...
        any(~isfinite(soundSquared),'all') || any(p<=0,'all') || ...
        any(cv<=0,'all') || any(pRho<=0,'all') || any(cp<=0,'all') || ...
        any(soundSquared<=0,'all')
    error('gasprop:UnstableTabularData', ...
        ['Tabular Gas spline has an invalid or unstable state inside a cell; ' ...
         'refine or restrict the single-phase table.']);
end
end
