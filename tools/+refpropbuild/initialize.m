function info=initialize(installPath)
%INITIALIZE Load REFPROP for offline table generation only.
% The supplied 64-bit prototype refers to a thunk DLL absent in some installs.
% Use the NIST 9.1.1 thunk kept beside this helper; leave the supplied files intact.
if nargin<1,installPath='D:\REFPROP';end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'validation','external','refprop'));
if ~libisloaded('refprop')
    prototype=@nativePrototype;
    loadlibrary(fullfile(installPath,'REFPRP64.DLL'),prototype,'alias','refprop');
    info=struct('FluidType','none','BasePath',[char(installPath) filesep], ...
        'FluidDir',['fluids' filesep],'nComp',0,'mixFlag',0,'z_mix',0);
    setappdata(0,'RefpropLoadedState',info);
end
info=getappdata(0,'RefpropLoadedState');
if ~isstruct(info) || ~isfield(info,'BasePath')
    error('refpropbuild:LibraryState','REFPROP is already loaded without the refpropm state; unload it before building.');
end
if ~strcmpi(regexprep(char(info.BasePath),'[\\/]+$',''),regexprep(char(installPath),'[\\/]+$',''))
    error('refpropbuild:InstallationConflict','Loaded REFPROP installation differs from InstallPath; unload it before building.');
end
end

function [methodinfo,structs,enuminfo,ThunkLibName]=nativePrototype()
[methodinfo,structs,enuminfo,ThunkLibName]=rp_proto64(fileparts(mfilename('fullpath')));
% Old MATLAB generated "long" names; NIST's 2016 thunk uses int32/cstring.
% Keep the user's pointer argument types and match only the thunk ABI name.
for i=1:numel(methodinfo.name)
    types=methodinfo.RHS{i};
    if isempty(types),types={};end
    types=regexprep(types,'^(double|long)Ptr$','voidPtr');
    types=regexprep(types,'^int8Ptr$','cstring');
    types=regexprep(types,'^long$','int32');
    methodinfo.thunkname{i}=['void' strjoin(types,'') 'Thunk'];
end
end
