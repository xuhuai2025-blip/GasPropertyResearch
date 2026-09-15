function outputs=buildRefpropTables(options)
%BUILDREFPROPTABLES Explicit offline REFPROP generation; never called by +gasprop.
%   addpath('tools'); outputs=buildRefpropTables;
%   buildRefpropTables(TemperatureCount=241,VolumeCount=641)
%   buildRefpropTables(Gases=["CO2","Air"]) % reference data still covers all gases
%   buildRefpropTables(SaveSource=true) % optional reproducibility intermediates
% Keep the supplied REFPROP files in validation/external/refprop. Main runtime
% uses generated MAT files and does not load REFPROP or add validation paths.
arguments
    options.InstallPath (1,1) string="D:\REFPROP"
    options.TemperatureCount (1,1) double=161
    options.VolumeCount (1,1) double=481
    options.PressureCount (1,1) double=161
    options.ReferenceCount (1,1) double=401
    options.Gases (1,:) string=["He","H2","N2","Ar","Air","CO2"]
    options.Compile (1,1) logical=true
    options.SaveSource (1,1) logical=false
    options.OutputDirectory (1,1) string=""
end
started=tic;
root=fileparts(fileparts(mfilename('fullpath')));
originalPath=path;
addpath(root);
originalValidation=gasprop.validation();
pathCleanup=onCleanup(@()restoreSession(originalValidation,originalPath)); %#ok<NASGU>
gasprop.validation(true);
gaspropcheck.refpropBuild('options',options);
installation=refpropbuild.initialize(options.InstallPath);
if strlength(options.OutputDirectory)==0,options.OutputDirectory=string(fullfile(root,'data'));end
if ~isfolder(options.OutputDirectory),mkdir(options.OutputDirectory);end
provenance=struct('GeneratedAtUTC',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')), ...
    'Generator','tools/buildRefpropTables.m','MATLABVersion',version, ...
    'REFPROPInstallation',installation.BasePath, ...
    'DLLSHA256',refpropbuild.sha256(fullfile(installation.BasePath,'REFPRP64.DLL')), ...
    'InterfaceSHA256',refpropbuild.sha256(fullfile(root,'validation','external','refprop','refpropm.m')), ...
    'ThunkSHA256',refpropbuild.sha256(fullfile(root,'tools','+refpropbuild','REFPRP64_thunk_pcwin64.dll')), ...
    'ThunkSource','https://trc.nist.gov/refprop/FAQ/MATLAB/9.1.1/REFPRP64_thunk_pcwin64.dll', ...
    'LegacyAPIDocumentation','https://refprop-docs.readthedocs.io/en/latest/DLL/legacy.html', ...
    'IdealCpMethod','THERM0dll Cp00 (J/mol/K) / REFPROP molar mass (kg/mol)', ...
    'ReferenceTransportMethod','TRNPRPdll(T,D=0); viscosity microPa s converted to Pa s');
fluids=refpropbuild.fluidConfig();
if any(~ismember(lower(options.Gases),lower(string({fluids.Key}))))
    error('refpropbuild:UnknownFluid','Gases must use configured keys: %s.',strjoin(string({fluids.Key}),', '));
end
referenceItems=cell(size(fluids));
infos=cell(size(fluids));
for i=1:numel(fluids)
    fluid=fluids(i);info=refpropbuild.fluidInfo(fluid,installation);
    referenceTemperature=max(300,info.MinimumTemperatureK);
    T=refpropbuild.uniqueGrid([logspace(log10(info.MinimumTemperatureK),log10(info.MaximumTemperatureK),options.ReferenceCount),referenceTemperature], ...
        [info.MinimumTemperatureK,info.MaximumTemperatureK]).';
    T([1 end])=[info.MinimumTemperatureK info.MaximumTemperatureK];
    [cp,mu,k]=refpropbuild.idealReference(T,info.REFPROPMolarMassKgPerMol);
    metadata=provenance;metadata.Fluid=info;
    referenceItems{i}=struct('TemperatureK',T,'IdealCpJPerKgK',cp, ...
        'ViscosityPaS',mu,'ThermalConductivityWPerMK',k, ...
        'ReferenceTemperatureK',referenceTemperature,'GasName',fluid.Name,'GasConstant',fluid.GasConstant, ...
        'MinimumTemperatureK',T(1),'MaximumTemperatureK',T(end),'PressureRangePa',[1e5 2e7], ...
        'CriticalTemperatureK',info.CriticalTemperatureK,'CriticalPressurePa',info.CriticalPressurePa, ...
        'SourceMetadata',metadata);
    fprintf('%s: reference T=[%g,%g] K, physical Tc=%g K\n',fluid.Name,T(1),T(end),info.CriticalTemperatureK);
    infos{i}=info;
end
rkReference=[referenceItems{:}];
rkFile=fullfile(options.OutputDirectory,'rk-reference.mat');
save(rkFile,'rkReference','-v7');
compiledFiles=strings(0,1);sourceFiles=strings(0,1);
for i=1:numel(fluids)
    fluid=fluids(i);
    if ~any(strcmpi(options.Gases,fluid.Key)),continue,end
    % Reference sampling leaves another fluid loaded; select this fluid again.
    refpropm('D','T',infos{i}.MinimumTemperatureK,'P',100,fluid.File);
    raw=refpropbuild.gasTable(fluid,infos{i},options,provenance);
    raw.SourceMetadata.GenerationWallTimeSeconds=toc(started);
    sourceName=string(fluid.OutputStem)+" REFPROP sampling (embedded provenance)";
    if options.SaveSource || ~options.Compile
        sourceFile=fullfile(options.OutputDirectory,[fluid.OutputStem '-refprop-source.mat']);
        save(sourceFile,'raw','-v7');
        sourceFiles(end+1)=sourceFile;
        sourceName=sourceFile;
    end
    if options.Compile
        addpath(root);
        gas=struct('Name',fluid.Name,'R',fluid.GasConstant,'MolarMass',fluid.MolarMass);
        tabularGas=gasprop.table.build(raw,gas,sourceName);
        compiledFile=fullfile(options.OutputDirectory,[fluid.OutputStem '-refprop.mat']);
        save(compiledFile,'tabularGas','-v7');
        compiledFiles(end+1)=compiledFile;
        fprintf('%s saved: %d x %d EOS, %d x %d transport\n',fluid.Key, ...
            size(raw.Compressibility),size(raw.TransportTable.ViscosityPaS));
    end
end
outputs=struct('ReferenceFile',rkFile,'SourceFiles',sourceFiles,'CompiledFiles',compiledFiles, ...
    'Compiled',options.Compile,'WallTimeSeconds',toc(started));
fprintf('%s\n',jsonencode(outputs));
end

function restoreSession(enabled,originalPath)
gasprop.validation(enabled);
path(originalPath);
end
