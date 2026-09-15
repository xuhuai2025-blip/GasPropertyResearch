function [p,pT,pRho]=pressure(table,R,T,rho)
%PRESSURE Mechanical EOS only, used in density inversion without caloric work.
v=1./rho;
[z,zT,zV]=gasprop.eos.internal.tabular.interpolate(table.CompressibilitySurface,v,T);
p=R.*T.*z./v;
pT=R.*(z+T.*zT)./v;
pRho=R.*T.*(z-v.*zV);
end
