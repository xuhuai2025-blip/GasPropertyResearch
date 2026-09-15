function entropy=deriveEntropy(T,v,z,R,referenceCv,zSurface,energySurface)
%DERIVEENTROPY Build s(v,T) from SAGE equation 18.38.
entropy=zeros(size(z));
entropy(:,end)=cumtrapz(T,referenceCv./T);
for temperatureIndex=1:numel(T)
    for volumeIndex=numel(v)-1:-1:1
        increment=intervalIntegral(zSurface,energySurface, ...
            T(temperatureIndex),v(volumeIndex),v(volumeIndex+1),R);
        entropy(temperatureIndex,volumeIndex)= ...
            entropy(temperatureIndex,volumeIndex+1)-increment;
    end
end
entropy=entropy+max(0,R-min(entropy(:)));
end

function increment=intervalIntegral(zSurface,energySurface,T,vLower,vUpper,R)
nodes=[-0.861136311594053,-0.339981043584856, ...
    0.339981043584856,0.861136311594053];
weights=[0.347854845137454,0.652145154862546, ...
    0.652145154862546,0.347854845137454];
halfWidth=(vUpper-vLower)/2;
v=(vUpper+vLower)/2+halfWidth*nodes;
z=gasprop.eos.internal.tabular.interpolate(zSurface,v,T+zeros(size(v)));
[~,~,energyV]=gasprop.eos.internal.tabular.interpolate( ...
    energySurface,v,T+zeros(size(v)));
increment=halfWidth*sum(weights.*(energyV+R.*T.*z./v)./T);
end
