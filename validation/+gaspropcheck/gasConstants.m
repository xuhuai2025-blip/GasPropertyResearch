function gasConstants(name,values)
%GASCONSTANTS Optional constructor checks for a custom gas.
validateattributes(string(name),{'string'},{'scalar','nonempty'});
validateattributes(values,{'numeric'},{'real','finite','positive','vector','numel',7});
if values(3)<=1
    error('gasprop:InvalidGamma','Specific heat ratio must be greater than one.');
end
end
