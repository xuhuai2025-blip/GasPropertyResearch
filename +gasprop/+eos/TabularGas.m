classdef TabularGas
    %TABULARGAS SAGE-style real gas from Z(v,T), epsilon(v,T), s(v,T).
    %   The compact EOS facade delegates loading, derived tables, spline
    %   interpolation, and safeguarded inversion to the tabular package.

    properties (SetAccess=private)
        R
        StateTable
    end

    methods
        function obj=TabularGas(gas,source)
            obj.R=gas.R;
            obj.StateTable=gasprop.eos.tabular.loadStateTable( ...
                source,gas);
        end

        function state=state(obj,T,rho)
            gasprop.eos.validateInputs(T,rho);
            state=gasprop.eos.tabular.state( ...
                obj.StateTable,obj.R,T,rho);
        end

        function rho=density(obj,T,p)
            gasprop.eos.validateInputs(T,p);
            rho=gasprop.eos.tabular.density( ...
                obj.StateTable,obj.R,T,p);
        end

        function value=metadata(obj)
            table=obj.StateTable;
            value=struct('Kind',table.Kind,'Source',table.Source, ...
                'SourceFile',table.SourceFile,'GasName',table.GasName, ...
                'GasConstant',table.GasConstant, ...
                'TemperatureRangeK',[table.TemperatureK(1),table.TemperatureK(end)], ...
                'DensityRangeKgPerM3',[table.MinimumDensityKgPerM3, ...
                    table.MaximumDensityKgPerM3], ...
                'ReferenceCpSource',table.ReferenceCpSource);
        end
    end
end
