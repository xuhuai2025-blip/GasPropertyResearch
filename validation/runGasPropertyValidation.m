function report=runGasPropertyValidation()
%RUNGASPROPERTYVALIDATION Optional, pure-property regression checks.
%   This file is deliberately outside +gasprop and is never invoked by the
%   package. Add the library root and this validation folder to path first.

if exist('gasprop.GasProperties','class')~=8
    error('gasprop:validation:LibraryNotOnPath', ...
        'Add the GasPropertyResearch root to the MATLAB path first.');
end
verifyModelCatalog();
idealError=verifyIdealGas();
rkError=verifyRKGas();
helium=verifyHelium();
tabular=verifyTabularGas();
report=struct('Passed',true, ...
    'IdealPressureRelativeError',idealError, ...
    'RKDensityRelativeError',rkError, ...
    'HeliumDensityRelativeError',helium.DensityRelativeError, ...
    'HeliumViscosityPaS',helium.ViscosityPaS, ...
    'TabularPressureRelativeError',tabular.PressureRelativeError, ...
    'TabularDensityRelativeError',tabular.DensityRelativeError, ...
    'TabularDerivativeRelativeError',tabular.DerivativeRelativeError);
fprintf(['GasPropertyResearch validation passed: ideal %.3g, RK density %.3g, ' ...
    'He density %.3g, tabular pressure %.3g.\n'], ...
    idealError,rkError,helium.DensityRelativeError, ...
    tabular.PressureRelativeError);
end

function verifyModelCatalog()
catalog=gasprop.modelCatalog();
expectedNames=["Ideal Gas","RK Gas","Tabular Gas", ...
    "NIST He-4 Helmholtz"];
expectedIds=["ideal","rk","tabular","helmholtz"];
require(isequal(string({catalog.Name}),expectedNames) && ...
    isequal(string({catalog.Id}),expectedIds), ...
    'ModelCatalog','The public four-model catalog changed unexpectedly.');
end

function relativeError=verifyIdealGas()
gas=gasprop.create("N2","ideal");
T=[250,300,500]; rho=[1.2,4.5,12.0];
state=gas.state(T,rho);
expected=rho.*gas.R.*T;
relativeError=max(abs(state.Pressure-expected)./expected);
require(relativeError<1e-12,'IdealPressure', ...
    'Ideal Gas pressure does not satisfy p = rho*R*T.');
recovered=gas.density(T,state.Pressure);
require(max(abs(recovered-rho)./rho)<1e-12,'IdealDensity', ...
    'Ideal Gas T-p density inversion did not recover rho.');
end

function relativeError=verifyRKGas()
gas=gasprop.create("N2","rk");
T=300; rho=5;
state=gas.state(T,rho);
recovered=gas.density(T,state.Pressure);
relativeError=abs(recovered-rho)/rho;
require(isfinite(state.Cp) && state.Cp>state.Cv && state.SoundSpeed>0, ...
    'RKState','RK Gas did not return a stable thermodynamic state.');
require(relativeError<1e-10,'RKDensity', ...
    'RK Gas T-p density inversion did not recover rho.');
end

function result=verifyHelium()
gas=gasprop.create("He","NIST He-4 Helmholtz","","nist-helium");
T=300; rho=2;
state=gas.state(T,rho);
recovered=gas.density(T,state.Pressure);
[mu,k,pr]=gas.transportAtDensity(T,rho);
relativeError=abs(recovered-rho)/rho;
require(all(isfinite([state.Pressure,state.Cv,state.Cp,state.SoundSpeed,mu,k,pr])), ...
    'HeliumFinite','NIST He-4 state or transport returned a non-finite value.');
require(state.Pressure>0 && state.Cv>0 && state.Cp>0 && ...
    state.SoundSpeed>0 && mu>0 && k>0 && pr>0, ...
    'HeliumPositive','NIST He-4 state or transport was not positive.');
require(relativeError<1e-10,'HeliumDensity', ...
    'NIST He-4 T-p density inversion did not recover rho.');
result=struct('DensityRelativeError',relativeError,'ViscosityPaS',mu);
end

function result=verifyTabularGas()
gas=gasprop.create("N2","ideal");
temperatureK=linspace(250,600,9).';
specificVolumeM3PerKg=logspace(log10(0.04),log10(2),13);
z=ones(numel(temperatureK),numel(specificVolumeM3PerKg));
table=buildTabularGasTable("N2",temperatureK, ...
    specificVolumeM3PerKg,z,struct('Source',"analytic ideal-N2 fixture"));
tabular=gas.withEquationOfState("tabular",table);
T=400; rho=3;
state=tabular.state(T,rho);
expected=rho*gas.R*T;
pressureRelativeError=abs(state.Pressure-expected)/expected;
recovered=tabular.density(T,expected);
densityRelativeError=abs(recovered-rho)/rho;
derivativeRelativeError=abs(state.PressureDensityDerivative-gas.R*T)/(gas.R*T);
require(pressureRelativeError<1e-10 && densityRelativeError<1e-10 && ...
    derivativeRelativeError<1e-10 && abs(state.Cv-gas.Cv)/gas.Cv<1e-10, ...
    'TabularIdeal','Ideal Z(v,T) table did not reproduce the ideal-gas contract.');

folder=string(tempname); mkdir(folder);
cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
matFile=fullfile(folder,"n2-table.mat");
prebuilt=buildTabularGasTable("N2",temperatureK, ...
    specificVolumeM3PerKg,z,struct('OutputFile',matFile));
fromMat=gas.withEquationOfState("tabular",matFile);
matState=fromMat.state(T,rho);
require(abs(matState.Pressure-state.Pressure)/state.Pressure<1e-12, ...
    'TabularMat','Prebuilt MAT table did not reproduce the in-memory state.');

textFile=fullfile(folder,"n2-z.txt");
writeSageTable(textFile,temperatureK,specificVolumeM3PerKg,z);
warningState=warning('query','gasprop:TabularCpFallback');
warning('off','gasprop:TabularCpFallback');
warningCleanup=onCleanup(@()warning(warningState)); %#ok<NASGU>
fromText=gas.withEquationOfState("tabular",textFile);
textState=fromText.state(T,rho);
require(abs(textState.Pressure-expected)/expected<1e-10, ...
    'TabularText','SAGE-style text table did not reproduce the ideal-gas state.');

[vGrid,tGrid]=meshgrid(specificVolumeM3PerKg,temperatureK);
nonidealZ=1+0.002*(tGrid/300).^2.*vGrid.^2;
nonideal=gas.withEquationOfState("tabular", ...
    buildTabularGasTable("N2",temperatureK,specificVolumeM3PerKg,nonidealZ));
nonidealState=nonideal.state(T,rho);
deltaT=1e-3;
finiteDifference=(nonideal.state(T+deltaT,rho).Pressure- ...
    nonideal.state(T-deltaT,rho).Pressure)/(2*deltaT);
nonidealDerivativeError=abs(finiteDifference- ...
    nonidealState.PressureTemperatureDerivative)/ ...
    abs(nonidealState.PressureTemperatureDerivative);
require(nonidealDerivativeError<1e-6,'TabularDerivative', ...
    'Tabular Gas pressure-temperature derivative failed a finite-difference check.');
result=struct('PressureRelativeError',pressureRelativeError, ...
    'DensityRelativeError',densityRelativeError, ...
    'DerivativeRelativeError',nonidealDerivativeError, ...
    'TableKind',prebuilt.Kind);
end

function writeSageTable(fileName,T,v,z)
fileId=fopen(fileName,'w');
if fileId<0
    error('gasprop:validation:TableWrite','Could not write %s.',fileName);
end
cleanup=onCleanup(@()fclose(fileId)); %#ok<NASGU>
fprintf(fileId,'T %c v',char(92));
fprintf(fileId,'\t%.17g',v);
fprintf(fileId,'\n');
for index=1:numel(T)
    fprintf(fileId,'%.17g',T(index));
    fprintf(fileId,'\t%.17g',z(index,:));
    fprintf(fileId,'\n');
end
end

function require(condition,name,message)
if ~condition
    error(['gasprop:validation:' char(name)],'%s',message);
end
end
