% get parameter from the struct 'model'
% input: model, the struct named 'model'
%      : temp, temperature, one of [0, 10, 25, 40, 55]
%      : paraName, a string, could be one of 
%      : hysteresis related: 'GParam', 'M0Para', 'MPara'
%      : 
% output: theParam, the parameter you want

function theParam = getParaECM(paraName, temp, model)
theFields = fieldnames(model.(['T',num2str(temp)]));
match = strcmpi(paraName,theFields);
if ~match
    error('Parameter "%s" does not exist in model', paraName)
end
fieldName = char(theFields(match));

theParam = model.(['T',num2str(temp)]).(fieldName);

end

