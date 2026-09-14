function validateInputs(T,x)
%VALIDATEINPUTS Positive, finite real thermodynamic inputs.
validateattributes(T,{'numeric'},{'real','finite','positive','nonempty'});
validateattributes(x,{'numeric'},{'real','finite','positive','nonempty'});
end
