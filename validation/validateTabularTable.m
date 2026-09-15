function report=validateTabularTable(table,options)
%VALIDATETABULARTABLE Explicit cell-interior stability/thermodynamic audit.
%   REPORT = VALIDATETABULARTABLE(COMPILED_TABLE) samples three positions in
%   each T/v cell. It never runs during construction, loading or state calls.
%   This checks interpolation consistency, not accuracy against REFPROP.
%   OPTIONS: Fractions (default 3-point interior set), Tolerance (0.05),
%   EnergyIdentityTolerance (0.001, relative to pressure terms),
%   ErrorOnFailure (false), IncludePadding (false), TemperatureChunkSize (32).
%   IncludePadding audits the complete EOS support grid. This audit evaluates
%   native cells directly and leaves the session validation switch unchanged,
%   so state-domain guards do not exclude interpolation support samples.
if nargin<2, options=struct(); end
fractions=option(options,'Fractions',[0.211324865405187,0.5,0.788675134594813]);
tolerance=option(options,'Tolerance',0.05);
energyTolerance=option(options,'EnergyIdentityTolerance',0.001);
chunkSize=option(options,'TemperatureChunkSize',32);
validateattributes(fractions,{'numeric'},{'real','finite','vector','>',0,'<',1});
validateattributes(tolerance,{'numeric'},{'real','finite','scalar','nonnegative'});
validateattributes(energyTolerance,{'numeric'},{'real','finite','scalar','nonnegative'});
validateattributes(chunkSize,{'numeric'},{'real','finite','scalar','integer','positive'});
required={'TemperatureK','SpecificVolumeM3PerKg','GasConstant', ...
    'CompressibilitySurface','InternalEnergySurface','EntropySurface'};
if ~isstruct(table) || ~isscalar(table) || ~all(isfield(table,required))
    error('gasprop:validation:CompiledTableRequired', ...
        'Pass a compiled table returned by gasprop.table.build or gasprop.table.load.');
end
tNodes=table.TemperatureK(:); vNodes=table.SpecificVolumeM3PerKg(:);
validateattributes(tNodes,{'numeric'},{'real','finite','positive','increasing','nonempty'});
validateattributes(vNodes,{'numeric'},{'real','finite','positive','increasing','nonempty'});
tSamples=sort(reshape(tNodes(1:end-1)+diff(tNodes).*fractions(:).',[],1));
vSamples=sort(reshape(vNodes(1:end-1)+diff(vNodes).*fractions(:).',[],1)).';
report=struct('Passed',true,'SampleFractions',fractions,'Tolerance',tolerance, ...
    'EnergyIdentityTolerance',energyTolerance, ...
    'SampleCount',0,'OutsideOperatingPressureCount',0, ...
    'InvalidStateCount',0,'MinimumPressurePa',Inf,'MinimumCvJPerKgK',Inf, ...
    'MinimumPressureDensityDerivative',Inf,'MinimumCpJPerKgK',Inf, ...
    'MinimumSoundSpeedSquared',Inf,'FirstInvalidState',struct(), ...
    'EnergyVolumeDerivative',emptyMismatch(), ...
    'EnergyIdentityResidual',emptyMismatch(), ...
    'EntropyTemperatureDerivative',emptyMismatch(), ...
    'EntropyVolumeDerivative',emptyMismatch(), ...
    'MismatchScale',"max(1, abs(actual), abs(expected)); quantities in SI", ...
    'EnergyIdentityScale',"max(abs(T * dp/dT), abs(p)); u_v = T * dp/dT - p", ...
    'Scope',"Cell interiors within the declared operating pressure envelope");
if option(options,'IncludePadding',false)
    report.Scope="All EOS support-grid cell interiors, including pressure padding";
end
for first=1:chunkSize:numel(tSamples)
    [v,T]=meshgrid(vSamples,tSamples(first:min(first+chunkSize-1,numel(tSamples))));
    [z,zT,zV]=gasprop.eos.internal.tabular.evaluateSurface(table.CompressibilitySurface,T,v);
    [~,uT,uV]=gasprop.eos.internal.tabular.evaluateSurface(table.InternalEnergySurface,T,v);
    [~,sT,sV]=gasprop.eos.internal.tabular.evaluateSurface(table.EntropySurface,T,v);
    p=table.GasConstant.*T.*z./v;
    included=true(size(p));
    if isfield(table,'PressureRangePa') && ~isempty(table.PressureRangePa) && ...
            ~option(options,'IncludePadding',false)
        included=p>=table.PressureRangePa(1) & p<=table.PressureRangePa(2);
    end
    report.OutsideOperatingPressureCount=report.OutsideOperatingPressureCount+nnz(~included);
    if ~any(included,'all'), continue, end
    T=T(included); v=v(included); z=z(included); zT=zT(included); zV=zV(included);
    p=p(included); uT=uT(included); uV=uV(included); sT=sT(included); sV=sV(included);
    stability=validateSurfaceState(T,v,z,zT,zV,uT,table.GasConstant);
    report.SampleCount=report.SampleCount+numel(T);
    report.InvalidStateCount=report.InvalidStateCount+stability.InvalidStateCount;
    for name=["MinimumPressurePa","MinimumCvJPerKgK", ...
            "MinimumPressureDensityDerivative","MinimumCpJPerKgK","MinimumSoundSpeedSquared"]
        report.(name)=min(report.(name),stability.(name));
    end
    if isempty(fieldnames(report.FirstInvalidState)) && stability.InvalidStateCount>0
        report.FirstInvalidState=stability.FirstInvalidState;
    end
    report.EnergyVolumeDerivative=updateMismatch(report.EnergyVolumeDerivative, ...
        uV,table.GasConstant.*T.^2.*zT./v,T,v);
    expected=table.GasConstant.*T.^2.*zT./v;
    report.EnergyIdentityResidual=updateMismatch(report.EnergyIdentityResidual, ...
        uV,expected,T,v,max(abs(p+expected),abs(p)));
    report.EntropyTemperatureDerivative=updateMismatch(report.EntropyTemperatureDerivative, ...
        sT,uT./T,T,v);
    report.EntropyVolumeDerivative=updateMismatch(report.EntropyVolumeDerivative, ...
        sV,(uV+p)./T,T,v);
end
for name=["EnergyVolumeDerivative","EnergyIdentityResidual"]
    peak=report.(name);
    if isfinite(peak.TemperatureKAtMaximum)
        T=peak.TemperatureKAtMaximum; v=peak.SpecificVolumeM3PerKgAtMaximum;
        [z,zT]=gasprop.eos.internal.tabular.evaluateSurface(table.CompressibilitySurface,T,v);
        peak.PressurePaAtMaximum=table.GasConstant*T*z/v;
        peak.TemperatureTimesPressureTemperatureDerivativePaAtMaximum= ...
            table.GasConstant*T*(z+T*zT)/v;
        peak.AbsoluteMismatchAtMaximum=abs(peak.ActualAtMaximum-peak.ExpectedAtMaximum);
        report.(name)=peak;
    end
end
report.LegacyRelativePassed=report.SampleCount>0 && report.InvalidStateCount==0 && ...
    report.EnergyVolumeDerivative.MaximumScaledMismatch<=tolerance && ...
    report.EntropyTemperatureDerivative.MaximumScaledMismatch<=tolerance && ...
    report.EntropyVolumeDerivative.MaximumScaledMismatch<=tolerance;
report.NaturalPressureResidual=report.EnergyIdentityResidual.MaximumScaledMismatch;
report.EnergyCriterionExplanation= ...
    "u_v = T*p_T-p subtracts two pressure terms. Relative error in u_v can diverge near cancellation; the acceptance residual uses max(abs(T*p_T),abs(p)). The former derivative-relative metric is retained unchanged for audit.";
report.Passed=report.SampleCount>0 && report.InvalidStateCount==0 && ...
    report.NaturalPressureResidual<=energyTolerance && ...
    report.EntropyTemperatureDerivative.MaximumScaledMismatch<=tolerance && ...
    report.EntropyVolumeDerivative.MaximumScaledMismatch<=tolerance;
fprintf(['Tabular audit: %d interior states, %d invalid; scaled mismatches ' ...
    'u_v %.4g, s_T %.4g, s_v %.4g; pressure residual %.4g; passed=%d.\n'], ...
    report.SampleCount,report.InvalidStateCount, ...
    report.EnergyVolumeDerivative.MaximumScaledMismatch, ...
    report.EntropyTemperatureDerivative.MaximumScaledMismatch, ...
    report.EntropyVolumeDerivative.MaximumScaledMismatch, ...
    report.NaturalPressureResidual,report.Passed);
if ~report.Passed && option(options,'ErrorOnFailure',false)
    error('gasprop:validation:TabularAuditFailed', ...
        ['Explicit table audit failed: %d invalid states; pressure-scaled energy ' ...
         'residual %.4g (limit %.4g); entropy [s_T,s_v] mismatches [%.4g,%.4g] (limit %.4g).'], ...
        report.InvalidStateCount,report.NaturalPressureResidual,energyTolerance, ...
        report.EntropyTemperatureDerivative.MaximumScaledMismatch, ...
        report.EntropyVolumeDerivative.MaximumScaledMismatch,tolerance);
end
end

function report=validateSurfaceState(T,v,z,zT,zV,cv,R)
%VALIDATESURFACESTATE Optional stability diagnostics for already sampled cells.
rho=1./v;
p=R.*T.*z./v;
pT=R.*(z+T.*zT)./v;
pRho=R.*T.*(z-v.*zV);
cp=cv+T.*pT.^2./(rho.^2.*pRho);
soundSquared=pRho+T.*pT.^2./(rho.^2.*cv);
invalid=~isfinite(p) | ~isfinite(cv) | ~isfinite(pRho) | ...
    ~isfinite(cp) | ~isfinite(soundSquared) | ...
    p<=0 | cv<=0 | pRho<=0 | cp<=0 | soundSquared<=0;
report=struct('InvalidStateCount',nnz(invalid), ...
    'MinimumPressurePa',min(p),'MinimumCvJPerKgK',min(cv), ...
    'MinimumPressureDensityDerivative',min(pRho), ...
    'MinimumCpJPerKgK',min(cp),'MinimumSoundSpeedSquared',min(soundSquared), ...
    'FirstInvalidState',struct());
if any(invalid)
    index=find(invalid,1);
    report.FirstInvalidState=struct('TemperatureK',T(index), ...
        'SpecificVolumeM3PerKg',v(index),'PressurePa',p(index), ...
        'CvJPerKgK',cv(index),'CpJPerKgK',cp(index), ...
        'PressureDensityDerivative',pRho(index), ...
        'SoundSpeedSquared',soundSquared(index));
end
end

function value=option(options,name,default)
value=default;
if isfield(options,name), value=options.(name); end
end

function result=emptyMismatch()
result=struct('MaximumScaledMismatch',0,'MaximumAbsoluteMismatch',0, ...
    'ActualAtMaximum',NaN,'ExpectedAtMaximum',NaN, ...
    'TemperatureKAtMaximum',NaN,'SpecificVolumeM3PerKgAtMaximum',NaN);
end

function result=updateMismatch(result,actual,expected,T,v,scale)
absolute=abs(actual-expected);
if nargin<6, scale=max(1,max(abs(actual),abs(expected))); end
scaled=absolute./scale;
scaled(~isfinite(scaled))=Inf;
[worst,index]=max(scaled);
result.MaximumAbsoluteMismatch=max(result.MaximumAbsoluteMismatch,max(absolute));
if worst>=result.MaximumScaledMismatch
    result.MaximumScaledMismatch=worst;
    result.ActualAtMaximum=actual(index); result.ExpectedAtMaximum=expected(index);
    result.TemperatureKAtMaximum=T(index); result.SpecificVolumeM3PerKgAtMaximum=v(index);
end
end
