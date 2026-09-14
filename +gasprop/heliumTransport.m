function [mu,k] = heliumTransport(T,rho)
%HELIUMTRANSPORT Helium viscosity [Pa s] and conductivity [W/(m K)].
% Viscosity: NIST TN 1334, pp. 4-5. Conductivity: Hands-Arp (1981),
% as implemented in CoolProp TransportRoutines.cpp (helium hardcoded).
% Background conductivity only; this implementation excludes T < 20 K.
% rho is kg/m^3. EOS and caloric properties are intentionally separate.
gasprop.validateHeliumTransportInputs(T,rho);
r=rho/1000; x=log(min(T,300));
low=exp(-.135311743./x+1.00347841+1.20654649*x ...
    -.149564551*x.^2+.0125208416*x.^3);
B=-47.5295259./x+87.6799309-42.0741589*x+8.33128289*x.^2-.589252385*x.^3;
C=547.309267./x-904.870586+431.404928*x-81.4504854*x.^2+5.37008433*x.^3;
D=-1684.39324./x+3331.08630-1632.19172*x+308.804413*x.^2-20.2936367*x.^3;
excess=low.*expm1(r.*B+r.^2.*C+r.^3.*D);
high=196*T.^.71938.*exp(12.451./T-295.67./T.^2-4.1249);
blend=min(max((T-100)/10,0),1);
mu=((1-blend).*low+blend.*high+excess)*1e-7;
s=3.739232544./T-26.20316969./T.^2+59.82252246./T.^3-49.26397634./T.^4;
k0=2.7870034e-3*T.^.7034007057.*exp(s);
t1=T.^(1/3); t2=T.^(2/3);
kr=(1.862970530e-4-7.275964435e-7*T-1.427549651e-4*t1+3.290833592e-5*t2).*rho ...
    +(-5.213335363e-8+4.492659933e-8*t1-5.924416513e-9*t2).*rho.^3 ...
    +(7.087321137e-6-6.013335678e-6*t1+8.067145814e-7*t2+3.995125013e-7./T) ...
    .*rho.^2.*log(max(rho,realmin)/68);
k=k0+kr;
end
