function tabularGas=buildTabularGasTable(gasName,temperatureK,specificVolumeM3PerKg,compressibility,options)
%BUILDTABULARGASTABLE Build a portable SAGE-style Tabular Gas state table.
%   This optional utility is outside the core package and never runs
%   automatically. COMPRESSIBILITY is temperature-by-volume.
%
%   OPTIONS may contain ReferenceCpJPerKgK, InternalEnergyJPerKg,
%   EntropyJPerKgK, Source, and OutputFile.

if nargin<5, options=struct(); end
if ~isstruct(options) || ~isscalar(options)
    error('gasprop:validation:InvalidOptions', ...
        'options must be a scalar struct.');
end
gas=gasprop.GasProperties.fromName(gasName);
raw=struct('TemperatureK',temperatureK, ...
    'SpecificVolumeM3PerKg',specificVolumeM3PerKg, ...
    'Compressibility',compressibility,'GasName',gas.Name, ...
    'MolarMass',gas.MolarMass,'GasConstant',gas.R, ...
    'ReferenceCpJPerKgK',gas.Cp);
names=["ReferenceCpJPerKgK","InternalEnergyJPerKg", ...
    "EntropyJPerKgK","Source"];
for name=names
    if isfield(options,char(name)), raw.(char(name))=options.(char(name)); end
end
eos=gasprop.eos.TabularGas(gas,raw);
tabularGas=eos.StateTable;
if isfield(options,'OutputFile') && strlength(string(options.OutputFile))>0
    outputFile=string(options.OutputFile);
    folder=fileparts(outputFile);
    if strlength(folder)>0 && ~isfolder(folder), mkdir(folder); end
    save(outputFile,'tabularGas');
end
end
