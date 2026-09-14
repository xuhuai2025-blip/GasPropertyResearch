function derivative=differentiate(values,grid,dimension)
%DIFFERENTIATE Cache node derivatives from the same cubic spline family.
grid=double(grid(:));
if numel(grid)<4 || any(diff(grid)<=0)
    error('gasprop:InvalidTabularData', ...
        'Derivative grid requires at least four strictly ascending nodes.');
end
if dimension==1
    if size(values,1)~=numel(grid)
        error('gasprop:InvalidTabularData', ...
            'Derivative grid length must match the first data dimension.');
    end
    derivative=zeros(size(values));
    for index=1:size(values,2)
        splineForm=spline(grid,values(:,index).');
        coefficients=splineForm.coefs;
        derivative(1:end-1,index)=coefficients(:,end-1);
        delta=grid(end)-grid(end-1);
        order=splineForm.order;
        powers=(order-1):-1:1;
        derivative(end,index)=sum(coefficients(end,1:end-1).* ...
            powers.*delta.^(powers-1));
    end
elseif dimension==2
    derivative=gasprop.eos.tabular.differentiate( ...
        values.',grid,1).';
else
    error('gasprop:InvalidTabularData', ...
        'Only two-dimensional tabular-gas data are supported.');
end
end
