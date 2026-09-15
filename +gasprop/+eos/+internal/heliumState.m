function s=heliumState(T,rho)
%HELIUMSTATE Normal-fluid He-4 Helmholtz EOS, SI mass-specific properties.
% NIST IR 8474 (2023), equations 2-4 and Table 2; reference state in Eq. 3.
% Thermodynamics only. Transport correlations have separate validity limits.
shape=size(T+rho); %#ok<ELARLOG> Intentional implicit expansion to a shared shape.
T=reshape(T+zeros(shape),1,[]);rho=reshape(rho+zeros(shape),1,[]);
R=8.3144598/.004002602;Tc=5.1953;rhoc=17383.7*.004002602;
delta=rho/rhoc;tau=Tc./T;
persistent n d t ell eta beta gam eps0
if isempty(n)
    n=[.015559018 3.0638932 -4.2420844 .054418088 -.18971904 .087856262 2.2833566 -.53331595 -.53296502 .99444915 -.30078896 -1.6432563 ...
       .8029102 .026838669 .04687678 -.14832766 .03016211 -.019986041 .14283514 .007418269 -.22989793 .79224829 -.049386338]';
    d=[4 1 1 2 2 3 1 1 3 2 2 1 2 1 2 1 1 3 2 2 3 2 2]';
    t=[1 .425 .63 .69 1.83 .575 .925 1.585 1.69 1.51 2.9 .8 1.26 3.51 2.785 1 4.22 .83 1.575 3.447 .73 1.634 6.13]';
    ell=[0 0 0 0 0 0 1 2 2 1 2 1]';
    eta=[1.5497 9.245 4.76323 6.3826 8.7023 .255 .3523 .1492 .05 .1668 42.2358]';
    beta=[.2471 .0983 .1556 2.6782 2.7077 .6621 .1775 .4821 .3069 .1758 1357.6577]';
    gam=[3.15 2.54505 1.2513 1.9416 .5984 2.2282 1.606 3.815 1.61958 .6407 1.076]';
    eps0=[.596 .3423 .761 .9747 .5868 .5627 2.5346 3.6763 4.5245 5.039 .959]';
end
v=n.*exp(d.*log(delta)+t.*log(tau));
ld=d./delta;lt=t./tau;ldd=-d./delta.^2;ltt=-t./tau.^2;
ix=7:12;e=ell(ix);de=delta.^e;
v(ix,:)=v(ix,:).*exp(-de);ld(ix,:)=ld(ix,:)-e.*de./delta;
ldd(ix,:)=ldd(ix,:)-e.*(e-1).*de./delta.^2;
ix=13:23;dd=delta-eps0;tt=tau-gam;
v(ix,:)=v(ix,:).*exp(-eta.*dd.^2-beta.*tt.^2);
ld(ix,:)=ld(ix,:)-2*eta.*dd;lt(ix,:)=lt(ix,:)-2*beta.*tt;
ldd(ix,:)=ldd(ix,:)-2*eta;ltt(ix,:)=ltt(ix,:)-2*beta;
ad=sum(v.*ld,1);at=sum(v.*lt,1);add=sum(v.*(ld.^2+ldd),1);
att=sum(v.*(lt.^2+ltt),1);adt=sum(v.*ld.*lt,1);
z=1+delta.*ad;p=rho.*R.*T.*z;
pT=rho.*R.*(z-delta.*tau.*adt);prho=R.*T.*(1+2*delta.*ad+delta.^2.*add);
u=R.*T.*(1.5+.4674522201550815*tau+tau.*at);
cv=R.*(1.5-tau.^2.*att);cp=cv+T.*pT.^2./(rho.^2.*prho);
c2=prho+T.*pT.^2./(rho.^2.*cv);
s=struct('Pressure',reshape(p,shape),'InternalEnergy',reshape(u,shape), ...
    'Enthalpy',reshape(u+p./rho,shape),'Cv',reshape(cv,shape),'Cp',reshape(cp,shape), ...
    'PressureTemperatureDerivative',reshape(pT,shape),'PressureDensityDerivative',reshape(prho,shape), ...
    'SoundSpeed',reshape(sqrt(c2),shape),'Z',reshape(z,shape),'R',R);
end
