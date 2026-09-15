function report=runRefpropValidation(options)
%RUNREFPROPVALIDATION Explicit independent-grid audit of all six gas tables.
%   No files are written unless OutputFolder is supplied. OPTIONS: Gases,
%   RefpropPath, RelativeTolerance (0.005), IncludeRK (true), DeepTableAudit
%   (true), OutputFolder, ReportName, TemperatureK, PressurePa. TableFile
%   requires exactly one selected gas. RK errors have no accuracy gate.
if nargin<1, options=struct(); end
folder=fileparts(mfilename('fullpath')); root=fileparts(folder);
oldPath=path;
addpath(root);
oldValidation=gasprop.validation();
cleanup=onCleanup(@()restoreConfiguration(oldValidation,oldPath)); %#ok<NASGU>
gasprop.validation(true);
addpath(fullfile(root,'tools'),fullfile(folder,'external','refprop'));
refpropbuild.initialize(option(options,'RefpropPath','D:\REFPROP'));
names=string(option(options,'Gases',{'He','H2','N2','Ar','Air','CO2'}));
fluids=refpropbuild.fluidConfig(); catalog=gasprop.gasCatalog();
if isfield(options,'TableFile') && numel(names)~=1
    error('gasprop:validation:SingleTable','TableFile requires exactly one gas in Gases.');
end
tolerance=option(options,'RelativeTolerance',0.005);
validateattributes(tolerance,{'numeric'},{'scalar','real','finite','positive'});
results=cell(size(names));rkResults=cell(size(names));
for index=1:numel(names)
    selected=find(strcmpi({fluids.Key},names(index)),1);
    if isempty(selected),error('gasprop:validation:Gas','Unknown gas %s.',names(index));end
    fluid=fluids(selected);
    base=gasprop.GasProperties.fromName(names(index));
    entry=catalog(string({catalog.Name})==base.Name);
    tableFile=option(options,'TableFile',fullfile(root,'data',entry.TabularFile));
    table=gasprop.table.load(tableFile,base);
    gas=gasprop.create(names(index),'tabular',table,'tabular');
    metadata=gas.equationOfStateMetadata();
    [tc,pc]=refpropm('TP','C',0,' ',0,fluid.File);pc=pc*1000;
    [T,P]=independentSamples(metadata,tc,pc,options);
    samples=referenceStates(T,P,fluid.File,gas.R);
    item=compareGas(gas,samples,metadata.PressureRangePa,tolerance,true);
    item.GasName=names(index);item.TableFile=string(tableFile);item.Metadata=metadata;
    item.TableSHA256=refpropbuild.sha256(tableFile);
    item.TableGrid=[numel(table.TemperatureK),numel(table.SpecificVolumeM3PerKg)];
    item.SampleCount=numel(T);
    item.Sampling="Independent T-p grid including domain edges and near-critical states";
    if option(options,'DeepTableAudit',true)
        item.TabularConsistency=validateTabularTable(table);
        item.Passed=item.Passed && item.TabularConsistency.Passed;
    end
    results{index}=item;
    fprintf('Tabular %s: %d samples; max rho/Cp/sound/mu/k %.5g/%.5g/%.5g/%.5g/%.5g; passed=%d.\n', ...
        names(index),numel(T),item.AtTemperaturePressure.Density.MaximumRelativeError, ...
        item.AtTemperaturePressure.Cp.MaximumRelativeError, ...
        item.AtTemperaturePressure.SoundSpeed.MaximumRelativeError, ...
        item.AtTemperaturePressure.Viscosity.MaximumRelativeError, ...
        item.AtTemperaturePressure.ThermalConductivity.MaximumRelativeError,item.Passed);
    if option(options,'IncludeRK',true)
        rk=gasprop.create(names(index),'rk');metadata=rk.equationOfStateMetadata();
        [rkT,rkP]=independentSamples(metadata,tc,pc,struct());
        reference=referenceStates(rkT,rkP,fluid.File,rk.R);
        rkItem=compareGas(rk,reference,metadata.PressureRangePa,tolerance,false);
        rkItem.GasName=names(index);rkItem.Metadata=metadata;
        rkItem.AccuracyGateApplied=false;rkResults{index}=rkItem;
        fprintf('RK %s: density max relative error %.5g, Cp %.5g, failed queries %d.\n', ...
            names(index),rkItem.AtTemperaturePressure.Density.MaximumRelativeError, ...
            rkItem.AtTemperaturePressure.Cp.MaximumRelativeError,rkItem.FailedQueryCount);
    end
end
report=struct('CreatedAtUTC',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss')), ...
    'Source',"User supplied refpropm interface / installed REFPROP", ...
    'MATLABVersion',version,'RelativeTolerance',tolerance, ...
    'Tabular',[results{:}],'Passed',all(cellfun(@(x)x.Passed,results)));
report.RKReferenceSHA256=refpropbuild.sha256(fullfile(root,'data','rk-reference.mat'));
if option(options,'IncludeRK',true),report.RK=[rkResults{:}];end
outputFolder=string(option(options,'OutputFolder',''));
if strlength(outputFolder)>0
    if ~isfolder(outputFolder),mkdir(outputFolder);end
    reportName=string(option(options,'ReportName','refprop-validation'));
    report.OutputMat=string(fullfile(outputFolder,reportName+".mat"));
    report.OutputJson=string(fullfile(outputFolder,reportName+".json"));
    save(report.OutputMat,'report');
    fid=fopen(report.OutputJson,'w','n','UTF-8');
    if fid<0,error('gasprop:validation:ReportWrite','Cannot write %s.',report.OutputJson);end
    fileCleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
end
fprintf('REFPROP audit completed: %d gases; passed=%d.\n',numel(names),report.Passed);
end

function samples=referenceStates(T,p,fluid,runtimeR)
names=["Density","InternalEnergy","Enthalpy","Entropy","Cv","Cp", ...
    "SoundSpeed","Z","Viscosity","ThermalConductivity", ...
    "PressureTemperatureDerivative","PressureDensityDerivative"];
samples=struct('TemperatureK',T,'Pressure',p);
for name=names, samples.(name)=zeros(size(T)); end
samples.Quality=zeros(size(T));
samples.REFPROPZ=zeros(size(T));
for index=1:numel(T)
    [d,u,h,s,cv,cp,a,z,mu,k,pt,drdp,q]= ...
        refpropm('DUHSOCAZVL#RQ','T',T(index),'P',p(index)/1000,char(fluid));
    values=[d,u,h,s,cv,cp,a,p(index)/(d*runtimeR*T(index)),mu,k,pt*1000,1000/drdp];
    if ~isfinite(q) || q<1 || ~all(isfinite(values))
        error('gasprop:validation:NonGasReference', ...
            '%s REFPROP returned a liquid/two-phase flag Q=%g at T=%g K, p=%g Pa.', ...
            fluid,q,T(index),p(index));
    end
    for j=1:numel(names), samples.(names(j))(index)=values(j); end
    samples.Quality(index)=q;
    samples.REFPROPZ(index)=z;
end
samples.ZConvention="p/(rho * runtime R * T); raw REFPROP Z retained separately";
end

function result=compareGas(gas,reference,pressureRange,tolerance,accuracyGate)
fields=["Pressure","Density","InternalEnergy","Enthalpy","Entropy", ...
    "Cv","Cp","SoundSpeed","Z","PressureTemperatureDerivative", ...
    "PressureDensityDerivative","Viscosity","ThermalConductivity"];
count=numel(reference.TemperatureK);
actual=struct(); sameDensity=struct();
for name=fields
    actual.(name)=NaN(count,1); sameDensity.(name)=NaN(count,1);
end
errors=strings(count,1); densityErrors=strings(count,1);
interior=reference.Pressure>pressureRange(1) & reference.Pressure<pressureRange(2);
for index=1:count
    T=reference.TemperatureK(index); p=reference.Pressure(index);
    try
        state=gas.stateTP(T,p); [mu,k]=gas.transport(T,p);
        state.Viscosity=mu; state.ThermalConductivity=k;
        for name=fields, actual.(name)(index)=state.(name); end
    catch exception
        errors(index)=string(exception.identifier)+": "+string(exception.message);
    end
    if interior(index)
        try
            state=gas.state(T,reference.Density(index));
            [mu,k]=gas.transportAtDensity(T,reference.Density(index));
            state.Viscosity=mu; state.ThermalConductivity=k;
            for name=fields, sameDensity.(name)(index)=state.(name); end
        catch exception
            densityErrors(index)=string(exception.identifier)+": "+string(exception.message);
        end
    end
end
result=struct('Reference',reference,'QueryErrors',errors, ...
    'SameDensityQueryErrors',densityErrors,'FailedQueryCount',nnz(strlength(errors)>0), ...
    'FailedSameDensityQueryCount',nnz(strlength(densityErrors)>0), ...
    'SameDensityIncluded',interior, ...
    'SameDensityScope',"Interior pressures only; exact pressure endpoints use public T,p queries", ...
    'AtTemperaturePressure',struct(),'AtTemperatureDensity',struct(), ...
    'EnergyReferenceConvention',"Raw u/h/s retained. Tabular energy gate uses max(abs(u),Cv*T) and max(abs(h),Cp*T); tables share the REFPROP reference convention.");
for name=fields
    result.AtTemperaturePressure.(name)=metric(actual.(name),reference.(name), ...
        reference.TemperatureK,reference.Pressure,true(count,1));
    result.AtTemperatureDensity.(name)=metric(sameDensity.(name),reference.(name), ...
        reference.TemperatureK,reference.Pressure,interior);
end
for name=["InternalEnergy","Enthalpy"]
    if name=="InternalEnergy",capacity=reference.Cv;else,capacity=reference.Cp;end
    scale=max(abs(reference.(name)),capacity.*reference.TemperatureK);
    result.AtTemperaturePressure.(name).MaximumThermalScaledError= ...
        max(abs(actual.(name)-reference.(name))./scale);
    result.AtTemperatureDensity.(name).MaximumThermalScaledError= ...
        max(abs(sameDensity.(name)(interior)-reference.(name)(interior))./scale(interior));
end
result.Passed=result.FailedQueryCount==0 && result.FailedSameDensityQueryCount==0;
result.AccuracyGateApplied=accuracyGate;
if accuracyGate
    for name=["Pressure","Density","Cv","Cp","SoundSpeed","Z", ...
            "PressureTemperatureDerivative","PressureDensityDerivative", ...
            "Viscosity","ThermalConductivity"]
        result.Passed=result.Passed && ...
            result.AtTemperaturePressure.(name).MaximumRelativeError<=tolerance && ...
            result.AtTemperatureDensity.(name).MaximumRelativeError<=tolerance;
    end
    for name=["InternalEnergy","Enthalpy"]
        result.Passed=result.Passed && ...
            result.AtTemperaturePressure.(name).MaximumThermalScaledError<=tolerance && ...
            result.AtTemperatureDensity.(name).MaximumThermalScaledError<=tolerance;
    end
end
end

function result=metric(actual,reference,T,p,included)
absolute=abs(actual-reference); relative=absolute./max(abs(reference),eps);
selected=find(included & isfinite(relative));
if isempty(selected)
    maximumRelative=Inf; maximumAbsolute=Inf; index=1;
else
    [maximumRelative,position]=max(relative(selected)); index=selected(position);
    maximumAbsolute=max(absolute(selected));
end
result=struct('Actual',actual,'Reference',reference,'AbsoluteError',absolute, ...
    'RelativeError',relative,'MaximumAbsoluteError',maximumAbsolute, ...
    'MaximumRelativeError',maximumRelative,'MaximumErrorSampleIndex',index, ...
    'TemperatureKAtMaximum',T(index),'PressurePaAtMaximum',p(index), ...
    'ActualAtMaximum',actual(index),'ReferenceAtMaximum',reference(index));
end

function [T,p]=independentSamples(metadata,tc,pc,options)
tr=metadata.TemperatureRangeK;pr=metadata.PressureRangePa;
ft=[0,.0063,.0217,.0571,.1037,.1731,.2813,.4197,.5731,.7319,.8837,1];
fp=[0,.0193,.0719,.1397,.2273,.3317,.4571,.5939,.7391,.8813,1];
t=unique([tr(1)+diff(tr)*ft,tr(1)+[.00037,.0037,.0217,.137,.4197,1.317], ...
    tc*[1.005,1.02,1.1]]);
t=t(t>=tr(1) & t<=tr(2));
p0=unique([exp(log(pr(1))+log(pr(2)/pr(1))*fp), ...
    pc*[.85,.95,.991,1.003,1.021,1.077,1.25,1.5],pr]);
p0(abs(p0-pr(1))<=64*eps(pr(1)))=pr(1);
p0(abs(p0-pr(2))<=64*eps(pr(2)))=pr(2);
p0=unique(p0);
p0=p0(p0>=pr(1) & p0<=pr(2));
t=option(options,'TemperatureK',t);p0=option(options,'PressurePa',p0);
validateGrid(t,tr,'TemperatureK');validateGrid(p0,pr,'PressurePa');
[p,T]=meshgrid(p0,t);T=T(:);p=p(:);
end

function validateGrid(grid,range,name)
validateattributes(grid,{'numeric'},{'real','finite','positive','vector','nonempty'});
if any(grid<range(1) | grid>range(2))
    error('gasprop:validation:QueryRange','%s grid is outside the declared operating range.',name);
end
end

function value=option(options,name,default)
value=default; if isfield(options,name), value=options.(name); end
end

function restoreConfiguration(enabled,originalPath)
gasprop.validation(enabled);
path(originalPath);
end
