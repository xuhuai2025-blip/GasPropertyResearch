function data=rkReference(gasName)
%RKREFERENCE Read bundled offline RK reference data without external software.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
file=fullfile(root,'data','rk-reference.mat');
data=[];
if ~isfile(file), return; end
persistent cachedFile cachedData cachedStamp
entry=dir(file); stamp=[entry.datenum,entry.bytes];
if isempty(cachedData) || ~strcmp(cachedFile,file) || ~isequal(cachedStamp,stamp)
    payload=builtin('load',file,'rkReference');
    cachedData=payload.rkReference;
    cachedFile=file;
    cachedStamp=stamp;
end
for i=1:numel(cachedData)
    if string(cachedData(i).GasName)==string(gasName)
        data=cachedData(i); return
    end
end
end
