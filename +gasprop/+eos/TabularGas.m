classdef TabularGas
    %TABULARGAS SAGE-style real gas from Z(v,T), epsilon(v,T), s(v,T).
    %   The compact EOS facade delegates loading, derived tables, spline
    %   interpolation, and safeguarded inversion to the tabular package.

    properties (SetAccess=private)
        R
        StateTable
    end
    properties (Access=private)
        ViscosityInterpolant
        ConductivityInterpolant
    end

    methods
        function obj=TabularGas(gas,source)
            obj.R=gas.R;
            obj.StateTable=gasprop.table.load( ...
                source,gas);
            if obj.StateTable.HasTransport
                tr=obj.StateTable.TransportTable;
                grids={tr.TemperatureK,log(tr.PressurePa(:))};
                obj.ViscosityInterpolant=griddedInterpolant(grids,tr.ViscosityPaS,'linear','none');
                obj.ConductivityInterpolant=griddedInterpolant(grids,tr.ThermalConductivityWPerMK,'linear','none');
            end
        end

        function state=state(obj,T,rho)
            if gasprop.validation(),gaspropcheck.eosTabular('stateInput',obj.StateTable,T,rho);end
            state=gasprop.eos.internal.tabular.state( ...
                obj.StateTable,obj.R,T,rho);
            if gasprop.validation(),gaspropcheck.eosTabular('state',obj.StateTable,state);end
        end

        function rho=density(obj,T,p)
            if gasprop.validation(),gaspropcheck.eosTabular('densityInput',obj.StateTable,T,p);end
            rho=gasprop.eos.internal.tabular.density( ...
                obj.StateTable,obj.R,T,p);
        end

        function value=metadata(obj)
            table=obj.StateTable;
            value=struct('Name',"tabular",'Description',"Tabular Gas", ...
                'Kind',table.Kind,'Source',table.Source, ...
                'SourceFile',table.SourceFile,'GasName',table.GasName, ...
                'GasConstant',table.GasConstant, ...
                'TemperatureRangeK',[table.TemperatureK(1),table.TemperatureK(end)], ...
                'DensityRangeKgPerM3',[table.MinimumDensityKgPerM3, ...
                    table.MaximumDensityKgPerM3], ...
                'ReferenceCpSource',table.ReferenceCpSource, ...
                'PressureRangePa',table.PressureRangePa,'HasTransport',table.HasTransport);
            if isfield(table,'SourceMetadata'), value.SourceMetadata=table.SourceMetadata; end
        end

        function [mu,k]=transport(obj,T,rho,p)
            if ~obj.StateTable.HasTransport
                error('gasprop:MissingTransportData','This table does not contain transport data.');
            end
            if nargin<4
                if gasprop.validation(),gaspropcheck.eosTabular('stateInput',obj.StateTable,T,rho);end
                p=gasprop.eos.internal.tabular.pressure(obj.StateTable,obj.R,T,rho);
            end
            if gasprop.validation(),gaspropcheck.transport('tabular',obj.StateTable,T,p);end
            shape=size(T+p);
            T=T+zeros(shape); p=p+zeros(shape);
            bounds=obj.StateTable.TransportTable.PressurePa([1,end]);
            % Remove only inversion roundoff at the declared inclusive endpoints.
            p=min(max(p,bounds(1)),bounds(2));
            mu=reshape(obj.ViscosityInterpolant(T(:),log(p(:))),shape);
            k=reshape(obj.ConductivityInterpolant(T(:),log(p(:))),shape);
            if gasprop.validation(),gaspropcheck.transport('state',mu,k);end
        end
    end
end
