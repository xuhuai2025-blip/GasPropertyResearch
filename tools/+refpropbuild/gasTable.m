function raw=gasTable(fluid,info,options,provenance)
%GASTABLE Sample gas-only T-v thermodynamics and a bounded T-p transport grid.
T=refpropbuild.temperatureGrid(info,options.TemperatureCount);
T([1 end])=[info.MinimumTemperatureK info.MaximumTemperatureK];
% Small volume padding includes the requested TP edges despite roundoff.
rhoEdges=zeros(numel(T),2);
for i=1:numel(T)
    rhoEdges(i,1)=refpropm('D','T',T(i),'P',100,fluid.File);
    rhoEdges(i,2)=refpropm('D','T',T(i),'P',20000,fluid.File);
end
v=logspace(log10(1/max(rhoEdges(:))/1.0001),log10(1/min(rhoEdges(:))*1.0001),options.VolumeCount);
% Extra critical-density nodes resolve sharp but single-phase derivatives.
if info.MinimumTemperatureK<1.1*info.CriticalTemperatureK
    vc=1/info.CriticalDensityKgPerM3;
    extra=vc*logspace(log10(0.5),log10(2),161);
    v=refpropbuild.uniqueGrid([v extra(extra>v(1) & extra<v(end))],[v(1),v(end)]);
end
rho=1./v;
z=[1 zeros(1,19)];
tmelt=zeros(size(T));pHighest=zeros(size(T));
for i=1:numel(T)
    pHighest(i)=refpropm('P','T',T(i),'D',max(rho),fluid.File)*1000;
    [~,~,tm,ierr,herr]=calllib('refprop','MELTPdll',pHighest(i)/1000,z,0,0,32*ones(255,1),255);
    if ierr~=0,error('refpropbuild:MeltingBoundary','%s',strtrim(char(herr(:)')));end
    tmelt(i)=tm;
    if gasprop.validation(), gaspropcheck.refpropBuild('padding',info,T(i),max(rho),pHighest(i),tm); end
end
fields={'Compressibility','InternalEnergyJPerKg','EntropyJPerKgK','CvJPerKgK', ...
    'PressureTemperatureDerivativePaPerK','PressureDensityDerivativePaM3PerKg'};
raw=struct('TemperatureK',T,'SpecificVolumeM3PerKg',v);
for f=fields,raw.(f{1})=zeros(numel(T),numel(v));end
for i=1:numel(T)
    for j=1:numel(v)
        [p,u,s,cv,pT,drdp,q]=refpropm('PUSO#RQ','T',T(i),'D',rho(j),fluid.File);
        if gasprop.validation(), gaspropcheck.refpropBuild('node',T(i),rho(j),p,u,s,cv,pT,drdp,q); end
        raw.Compressibility(i,j)=p*1000/(rho(j)*fluid.GasConstant*T(i));
        raw.InternalEnergyJPerKg(i,j)=u;
        raw.EntropyJPerKgK(i,j)=s;
        raw.CvJPerKgK(i,j)=cv;
        raw.PressureTemperatureDerivativePaPerK(i,j)=pT*1000;
        raw.PressureDensityDerivativePaM3PerKg(i,j)=1000/drdp;
    end
    if mod(i,25)==0,fprintf('%s EOS: %d/%d temperature rows\n',fluid.Key,i,numel(T));end
end
P=logspace(5,log10(2e7),options.PressureCount);
if info.MinimumTemperatureK<1.1*info.CriticalTemperatureK
    pc=info.CriticalPressurePa;
    offset=logspace(-4,log10(2),161);
    extra=[pc*(1-offset),pc,pc*(1+offset)];
    P=refpropbuild.uniqueGrid([P extra(extra>1e5 & extra<2e7)],[1e5,2e7]);
end
P([1 end])=[1e5 2e7];
% Refuse silent extrapolation from the installed transport models.
if gasprop.validation(), gaspropcheck.refpropBuild('transportPressure',info,P(end)); end
transportT=T;
if info.MinimumTemperatureK<1.1*info.CriticalTemperatureK
    % Transport uses bilinear interpolation; refine both axes independently
    % of the Hermite EOS grid to resolve near-critical enhancement peaks.
    fraction=(1:fluid.TransportRefinement-1)/fluid.TransportRefinement;
    extraT=T(1:end-1)+diff(T)*fraction;
    logP=log(P);
    extraP=exp(logP(1:end-1).'+diff(logP).'*fraction);
    transportT=refpropbuild.uniqueGrid([T;extraT(:)],[T(1),T(end)]);
    P=refpropbuild.uniqueGrid([P,extraP(:).'],[1e5,2e7]);
end
mu=zeros(numel(transportT),numel(P));k=mu;
for i=1:numel(transportT)
    for j=1:numel(P)
        [mu(i,j),k(i,j),rhoTP,q]=refpropm('VLDQ','T',transportT(i),'P',P(j)/1000,fluid.File);
        if gasprop.validation(), gaspropcheck.refpropBuild('transportNode',info,transportT(i),P(j),mu(i,j),k(i,j),rhoTP,q); end
    end
    if mod(i,100)==0
        fprintf('%s transport: %d/%d temperature rows\n',fluid.Key,i,numel(transportT));
    end
end
raw.TransportTable=struct('TemperatureK',transportT,'PressurePa',P, ...
    'ViscosityPaS',mu,'ThermalConductivityWPerMK',k);
[cp0,~,~]=refpropbuild.idealReference(T,info.REFPROPMolarMassKgPerMol);
raw.ReferenceCpJPerKgK=cp0;
raw.GasName=fluid.Name;raw.GasConstant=fluid.GasConstant;raw.MolarMass=fluid.MolarMass;
raw.PressureRangePa=[1e5 2e7];
raw.Source=sprintf('REFPROP %g via user-supplied refpropm.m; %s',info.Version,fluid.File);
raw.SourceMetadata=provenance;
raw.SourceMetadata.Fluid=info;
raw.SourceMetadata.RequestedTemperatureRangeK=[10 800];
raw.SourceMetadata.RuntimeTemperatureRangeK=[T(1) T(end)];
raw.SourceMetadata.RuntimePressureRangePa=raw.PressureRangePa;
raw.SourceMetadata.EOSGridSize=[numel(T),numel(v)];
raw.SourceMetadata.TransportGridSize=size(mu);
raw.SourceMetadata.TransportRefinement=fluid.TransportRefinement;
raw.SourceMetadata.GridMethod='Log T and v; near-critical log(T-Tc), critical-density/pressure refinement; 64-ulp roundoff duplicates merged before interpolation.';
raw.SourceMetadata.SpecificVolumeRangeM3PerKg=[v(1) v(end)];
raw.SourceMetadata.MaximumPaddingPressurePa=max(pHighest);
raw.SourceMetadata.MaximumPressureByTemperaturePa=pHighest;
raw.SourceMetadata.MeltingTemperatureAtPaddingPressureK=tmelt;
raw.SourceMetadata.MeltingCheck='MELTP(p_max(T)); returned melting temperatures are inside the source model range and below T. No high-temperature MELTT extrapolation.';
raw.SourceMetadata.TransportCoordinates='T-p grid; no transport values evaluated at high-pressure T-v padding nodes.';
raw.SourceMetadata.ZDefinition='p/(rho*GasConstant*T); intentionally uses runtime R, not the REFPROP Z output.';
raw.SourceMetadata.DerivativeUnits='refpropm # * 1000 = dp/dT [Pa/K]; 1000 / refpropm R = dp/drho [Pa m3/kg]';
end
