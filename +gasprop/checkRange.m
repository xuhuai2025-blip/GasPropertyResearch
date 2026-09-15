function checkRange(value,bounds,identifier,label,tolerance)
%CHECKRANGE Always-on core range guard; no validation-folder dependency.
if nargin<5, tolerance=0; end
if isempty(bounds), bounds=[realmin Inf]; end
if ~isreal(value) || any(~isfinite(value(:)) | ...
        value(:)<bounds(1)*(1-tolerance) | value(:)>bounds(2)*(1+tolerance))
    error(identifier,'%s requires %.12g <= value <= %.12g.',label,bounds(1),bounds(2));
end
end
