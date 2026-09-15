function T=temperatureGrid(info,count)
%TEMPERATUREGRID Log-T baseline plus near-critical supercritical refinement.
lo=info.MinimumTemperatureK;hi=info.MaximumTemperatureK;
T=logspace(log10(lo),log10(hi),count);
T([1 end])=[lo hi];
if lo<1.1*info.CriticalTemperatureK
    tc=info.CriticalTemperatureK;
    extra=tc+logspace(log10(lo-tc),log10(hi-tc),count);
    T=[T extra(extra>lo & extra<hi)];
end
T=refpropbuild.uniqueGrid([T,max(300,lo)],[lo,hi]).';
end
