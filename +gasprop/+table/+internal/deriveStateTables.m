function [energy,entropy,energyVolumeDerivative]=deriveStateTables( ...
        T,v,z,R,cp0,zT,zSurface)
%DERIVESTATETABLES Apply SAGE 18.1.4 equations 18.36-18.38 once at load.
if nargin<6 || isempty(zT)
    zT=gasprop.table.internal.differentiate(z,T,1);
end
if nargin<7 || isempty(zSurface)
    zSurface=gasprop.table.internal.buildSurface(T,v,z,zT, ...
        gasprop.table.internal.differentiate(z,v,2));
end
energyVolumeDerivative=R.*(T.^2).*zT./v;
energy=zeros(size(z));
energy(:,end)=cumtrapz(T,cp0-R);
entropy=zeros(size(z));
entropy(:,end)=cumtrapz(T,(cp0-R)./T);
for temperatureIndex=1:numel(T)
    for volumeIndex=numel(v)-1:-1:1
        [energyIncrement,entropyIncrement]=intervalIntegrals( ...
            zSurface,T(temperatureIndex),v(volumeIndex),v(volumeIndex+1),R);
        energy(temperatureIndex,volumeIndex)= ...
            energy(temperatureIndex,volumeIndex+1)-energyIncrement;
        entropy(temperatureIndex,volumeIndex)= ...
            entropy(temperatureIndex,volumeIndex+1)-entropyIncrement;
    end
end
energyFloor=(cp0(1)-R)*T(1);
energy=energy+max(0,energyFloor-min(energy(:)));
entropy=entropy+max(0,R-min(entropy(:)));
end

function [energyIncrement,entropyIncrement]=intervalIntegrals(surface,T,vLower,vUpper,R)
nodes=[-0.861136311594053,-0.339981043584856, ...
    0.339981043584856,0.861136311594053];
weights=[0.347854845137454,0.652145154862546, ...
    0.652145154862546,0.347854845137454];
halfWidth=(vUpper-vLower)/2;
v=(vUpper+vLower)/2+halfWidth*nodes;
[z,zT]=gasprop.eos.internal.tabular.interpolate( ...
    surface,v,T+zeros(size(v)));
energyIncrement=halfWidth*sum(weights.*R.*T.^2.*zT./v);
entropyIncrement=halfWidth*sum(weights.*R.*(z+T.*zT)./v);
end
