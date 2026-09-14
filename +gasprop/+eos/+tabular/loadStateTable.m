function table=loadStateTable(source,gas)
%LOADSTATETABLE Load a SAGE-format text table or a prebuilt MAT state table.
if isstruct(source)
    raw=source;
    sourceName="in-memory table";
elseif isstring(source) || ischar(source)
    sourceName=string(source);
    if endsWith(lower(sourceName),".mat")
        payload=load(sourceName);
        if isfield(payload,'tabularGas')
            raw=payload.tabularGas;
        else
            raw=payload;
        end
    else
        raw=readSageCompressibility(sourceName);
    end
else
    error('gasprop:InvalidTabularData', ...
        'Tabular Gas source must be a table struct or a .mat/.txt file.');
end
table=gasprop.eos.tabular.normalizeStateTable( ...
    raw,gas,sourceName);
end

function raw=readSageCompressibility(fileName)
% SAGE 28.4: first row is ascending v; first column is ascending T.
values=readmatrix(fileName,'FileType','text');
if size(values,1)<5 || size(values,2)<5
    error('gasprop:InvalidTabularData', ...
        'SAGE Z(v,T) text table needs one header row/column and at least 4-by-4 data.');
end
raw=struct('TemperatureK',values(2:end,1), ...
    'SpecificVolumeM3PerKg',values(1,2:end), ...
    'Compressibility',values(2:end,2:end));
end
