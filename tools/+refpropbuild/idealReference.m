function [cp,mu,k]=idealReference(T,molarMass)
%IDEALREFERENCE Cp0 is exact ideal gas; transport is the zero-density limit.
% REFPROP's THERM0 returns J/(mol K), TRNPRP returns microPa s and W/(m K).
z=[1 zeros(1,19)];
cp=zeros(size(T));mu=cp;k=cp;
for i=1:numel(T)
    [~,~,~,~,~,~,~,~,cpMol]=calllib('refprop','THERM0dll',T(i),1e-8,z,0,0,0,0,0,0,0,0,0);
    [~,~,~,eta,k(i),ierr,herr]=calllib('refprop','TRNPRPdll',T(i),0,z,0,0,0,32*ones(255,1),255);
    if ierr~=0
        error('refpropbuild:ReferenceTransport','T=%g K: %s',T(i),strtrim(char(herr(:)')));
    end
    cp(i)=cpMol/molarMass;
    mu(i)=eta*1e-6;
end
if gasprop.validation(), gaspropcheck.refpropBuild('reference',cp(:),mu(:),k(:)); end
end
