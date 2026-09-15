function tableSurface(stage,values,first,second)
%TABLESURFACE Optional spline-grid and Hermite-node shape checks.
switch stage
    case 'surface'
        if ~isequal(size(values),size(first),size(second))
            error('gasprop:InvalidTabularData', ...
                'A tabular surface and its node gradients must have the same shape.');
        end
    case 'derivative'
        grid=first; dimension=second;
        if numel(grid)<4 || any(diff(grid)<=0)
            error('gasprop:InvalidTabularData', ...
                'Derivative grid requires at least four strictly ascending nodes.');
        end
        if ~isscalar(dimension) || ~any(dimension==[1 2])
            error('gasprop:InvalidTabularData', ...
                'Only two-dimensional tabular-gas data are supported.');
        end
        if size(values,dimension)~=numel(grid)
            error('gasprop:InvalidTabularData', ...
                'Derivative grid length must match the selected data dimension.');
        end
end
end
