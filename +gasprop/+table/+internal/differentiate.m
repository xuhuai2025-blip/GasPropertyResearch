function derivative=differentiate(values,grid,dimension)
%DIFFERENTIATE Cache node derivatives from the same cubic spline family.
grid=double(grid(:));
if gasprop.validation(), gaspropcheck.tableSurface('derivative',values,grid,dimension); end
if dimension==1
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
    derivative=gasprop.table.internal.differentiate( ...
        values.',grid,1).';
end
end
