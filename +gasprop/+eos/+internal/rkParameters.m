function parameters=rkParameters(gasName)
%RKPARAMETERS Published Sage c1/c2 fits; critical points may be overridden by
% REFPROP reference data. SageStlxHyperlinked(1).pdf, printed pp.323-324.
source="Sage manual chapter 28, printed pp.323-324 (PDF pp.327-328)";
pseudoPure=false;
switch string(gasName)
    case "Helium"
        Tc=5.195300013635951; pc=228322.7892147868; c1=.5338; c2=.05886;
    case "Hydrogen"
        Tc=33.1443326883113; pc=1296357.6060553084; c1=.4442; c2=.08207;
    case "Nitrogen"
        Tc=126.192; pc=3395800; c1=.4456; c2=.1102;
    case "Air"
        % Default pseudo-pure critical point from supplied REFPROP AIR.PPF
        % header. The actual REFPROP query result takes precedence when loaded.
        Tc=132.5306; pc=3786000; c1=.4372; c2=.1063; pseudoPure=true;
    case "Argon"
        Tc=150.687; pc=4863000; c1=.4274802335403414; c2=.08664034996495772;
        source="Standard RK critical-point constants; no Sage argon fit supplied";
    case "CarbonDioxide"
        Tc=304.1282; pc=7377300; c1=.4274802335403414; c2=.08664034996495772;
        source="Standard RK critical-point constants; no Sage CO2 fit supplied";
    otherwise
        error('gasprop:EOSGasCombination','RK supports He, H2, N2, Ar, Air and CO2.');
end
parameters=struct('CriticalTemperatureK',Tc,'CriticalPressurePa',pc, ...
    'C1',c1,'C2',c2,'Source',source,'IsPseudoPure',pseudoPure);
end
