function gas=create(gasName,equationOfState,tabularData,transportModel)
%CREATE Construct a reusable gasprop.GasProperties instance.
%   GAS = GASPROP.CREATE("He") creates the default He model.
%   GAS = GASPROP.CREATE("N2","rk") selects RK Gas.
%   GAS = GASPROP.CREATE("N2","tabular",TABLE) selects Tabular Gas.
%   Selecting Tabular Gas uses the selected gas's bundled REFPROP table.
%   RK defaults to bundled zero-density reference transport when available.
%   The optional fourth argument explicitly overrides transport selection.

if nargin<1 || strlength(string(gasName))==0
    error('gasprop:MissingGasName','A supported gas name is required.');
end
gas=gasprop.GasProperties.fromName(gasName);
if nargin>=2 && ~isempty(equationOfState)
    if nargin<3, tabularData=[]; end
    gas=gas.withEquationOfState(equationOfState,tabularData);
end
gas=gas.withDefaultTransport();
if nargin>=4 && ~isempty(transportModel)
    gas=gas.withTransportModel(transportModel);
end
end
