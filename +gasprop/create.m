function gas=create(gasName,equationOfState,tabularData,transportModel)
%CREATE Construct a reusable gasprop.GasProperties instance.
%   GAS = GASPROP.CREATE("He") creates the default He model.
%   GAS = GASPROP.CREATE("N2","rk") selects RK Gas.
%   GAS = GASPROP.CREATE("N2","tabular",TABLE) selects Tabular Gas.
%   The optional fourth argument selects sutherland or nist-helium transport.

if nargin<1 || strlength(string(gasName))==0
    error('gasprop:MissingGasName','A supported gas name is required.');
end
gas=gasprop.GasProperties.fromName(gasName);
if nargin>=2 && ~isempty(equationOfState)
    if nargin<3, tabularData=[]; end
    gas=gas.withEquationOfState(equationOfState,tabularData);
end
if nargin>=4 && ~isempty(transportModel)
    gas=gas.withTransportModel(transportModel);
end
end
