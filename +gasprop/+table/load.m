function table=load(source,gas)
%LOAD Read compiled data unchanged, or compile a raw SAGE T-v source once.
if gasprop.validation(), gaspropcheck.tableLoad('source',source); end
if isstruct(source)
    raw=source;
    sourceName="in-memory table";
else
    sourceName=string(source);
    if endsWith(lower(sourceName),".mat")
        payload=builtin('load',char(sourceName));
        if isfield(payload,'tabularGas')
            raw=payload.tabularGas;
        elseif isfield(payload,'raw')
            raw=payload.raw;
        else
            raw=payload;
        end
    else
        values=readmatrix(sourceName,'FileType','text');
        if gasprop.validation(), gaspropcheck.tableLoad('text',values); end
        raw=struct('TemperatureK',values(2:end,1), ...
            'SpecificVolumeM3PerKg',values(1,2:end), ...
            'Compressibility',values(2:end,2:end));
    end
end
if isstruct(raw) && isscalar(raw) && isfield(raw,'SchemaVersion')
    table=raw;
    if ~isfield(table,'PressureRangePa'), table.PressureRangePa=[]; end
    if ~isfield(table,'HasTransport'), table.HasTransport=false; end
    if gasprop.validation(), gaspropcheck.tableLoad('compiled',table,gas); end
else
    table=gasprop.table.build(raw,gas,sourceName);
end
end