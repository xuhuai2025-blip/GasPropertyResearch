function [value,temperatureDerivative,volumeDerivative]=evaluateSurface(surface,T,v)
%EVALUATESURFACE Evaluate a cached bicubic surface and its exact derivatives.
shape=size(T+v); %#ok<ELARLOG> Deliberate scalar expansion to common shape.
T=reshape(T+zeros(shape),[],1); v=reshape(v+zeros(shape),[],1);
tGrid=surface.TemperatureK; vGrid=surface.SpecificVolumeM3PerKg;
tCell=cellIndex(T,tGrid); vCell=cellIndex(v,vGrid);
tLeft=tGrid(tCell); tLeft=tLeft(:);
vLeft=vGrid(vCell); vLeft=vLeft(:);
deltaT=tGrid(tCell+1)-tLeft; deltaT=deltaT(:);
deltaV=vGrid(vCell+1); deltaV=deltaV(:)-vLeft;
eta=(T-tLeft)./deltaT; xi=(v-vLeft)./deltaV;
[h0x,h1x,g0x,g1x,dh0x,dh1x,dg0x,dg1x]=basis(xi,deltaV);
[h0t,h1t,g0t,g1t,dh0t,dh1t,dg0t,dg1t]=basis(eta,deltaT);
size2=[numel(tGrid),numel(vGrid)];
index00=sub2ind(size2,tCell,vCell);
index01=sub2ind(size2,tCell,vCell+1);
index10=sub2ind(size2,tCell+1,vCell);
index11=sub2ind(size2,tCell+1,vCell+1);
[f00,f01,f10,f11]=corners(surface.Value);
[fv00,fv01,fv10,fv11]=corners(surface.VolumeDerivative);
[ft00,ft01,ft10,ft11]=corners(surface.TemperatureDerivative);
[fvt00,fvt01,fvt10,fvt11]=corners(surface.CrossDerivative);
value=combine(h0x,h1x,g0x,g1x,h0t,h1t,g0t,g1t, ...
    f00,f01,f10,f11,fv00,fv01,fv10,fv11,ft00,ft01,ft10,ft11, ...
    fvt00,fvt01,fvt10,fvt11,deltaV,deltaT);
if nargout>2
volumeDerivative=combine(dh0x,dh1x,dg0x,dg1x,h0t,h1t,g0t,g1t, ...
    f00,f01,f10,f11,fv00,fv01,fv10,fv11,ft00,ft01,ft10,ft11, ...
    fvt00,fvt01,fvt10,fvt11,deltaV,deltaT);
volumeDerivative=reshape(volumeDerivative,shape);
end
if nargout>1
temperatureDerivative=combine(h0x,h1x,g0x,g1x,dh0t,dh1t,dg0t,dg1t, ...
    f00,f01,f10,f11,fv00,fv01,fv10,fv11,ft00,ft01,ft10,ft11, ...
    fvt00,fvt01,fvt10,fvt11,deltaV,deltaT);
temperatureDerivative=reshape(temperatureDerivative,shape);
end
value=reshape(value,shape);

    function [a00,a01,a10,a11]=corners(values)
        a00=values(index00); a01=values(index01);
        a10=values(index10); a11=values(index11);
        a00=a00(:); a01=a01(:); a10=a10(:); a11=a11(:);
    end
end

function index=cellIndex(query,grid)
index=sum(query>=grid(:).',2);
index=min(index,numel(grid)-1);
end

function [h0,h1,g0,g1,dh0,dh1,dg0,dg1]=basis(q,scale)
h0=2*q.^3-3*q.^2+1; h1=-2*q.^3+3*q.^2;
g0=q.^3-2*q.^2+q; g1=q.^3-q.^2;
dh0=(6*q.^2-6*q)./scale; dh1=(-6*q.^2+6*q)./scale;
dg0=(3*q.^2-4*q+1)./scale; dg1=(3*q.^2-2*q)./scale;
end

function value=combine(h0x,h1x,g0x,g1x,h0t,h1t,g0t,g1t, ...
        f00,f01,f10,f11,fv00,fv01,fv10,fv11,ft00,ft01,ft10,ft11, ...
        fvt00,fvt01,fvt10,fvt11,deltaV,deltaT)
value=h0x.*h0t.*f00+h1x.*h0t.*f01+h0x.*h1t.*f10+h1x.*h1t.*f11+ ...
    deltaV.*(g0x.*h0t.*fv00+g1x.*h0t.*fv01+g0x.*h1t.*fv10+g1x.*h1t.*fv11)+ ...
    deltaT.*(h0x.*g0t.*ft00+h1x.*g0t.*ft01+h0x.*g1t.*ft10+h1x.*g1t.*ft11)+ ...
    deltaV.*deltaT.*(g0x.*g0t.*fvt00+g1x.*g0t.*fvt01+ ...
        g0x.*g1t.*fvt10+g1x.*g1t.*fvt11);
end
