function eosInputs(varargin)
%EOSINPUTS Optional positive, finite, real thermodynamic input checks.
for index=1:nargin
    validateattributes(varargin{index},{'numeric'}, ...
        {'real','finite','positive','nonempty'});
end
end
