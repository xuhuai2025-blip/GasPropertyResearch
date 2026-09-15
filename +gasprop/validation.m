function enabled=validation(setting)
%VALIDATION Session-wide optional checks; disabled in a fresh MATLAB session.
%   gasprop.validation(true) enables checks in validation/+gaspropcheck.
%   gasprop.validation(false) disables them; gasprop.validation() reads the flag.
persistent active
if isempty(active), active=false; end
if nargin
    if ~islogical(setting) || ~isscalar(setting)
        error('gasprop:ValidationSetting','Validation setting must be true or false.');
    end
    if setting
        root=fileparts(fileparts(mfilename('fullpath')));
        addpath(fullfile(root,'validation'));
    end
    active=setting;
end
enabled=active;
end
