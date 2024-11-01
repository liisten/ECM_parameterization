%
%  input: struct, a struct like 'data', which contains dynamic test data.
%       : variable, time,current,voltage,capacity,soc and other...
%       : start,  2~8, 2 for 20%, 8 for 80%
%       : operationType: '模拟工况' '恒流充电' '恒功率充电' '高搁置' '低搁置'
%       : 高低指SOC高低，对充电，低搁置是充电前的搁置，高搁置是充电后的搁置
% output: theData, the data you want.
%       : index, 第start个开始的索引
%       : indexList, 第start个开始的连续索引

function theData = getDynData(struct,variable,start,operationType)
theFields = fieldnames(struct);

% check if the required variable exists
if ~ismember(theFields,variable)
    error('The required data does not exist: %s.',variable);
end

if ismember(operationType,{'模拟工况', '恒流充电', '恒功率充电'})
    % 找到连续样本中的第一个样本索引
    indices = find([false, diff(strcmp({struct.type}, operationType)) == 1]);
    indices = indices(:);
    index = indices(9-start);

    % 第start个开始的连续索引数组
    indexList = [];
    for k = index:length(struct)
        if strcmp(struct(k).type, operationType)
            indexList = [indexList;k];
        else
            break;
        end
    end

    theData = [];
    for k = 1:length(indexList)
        x = getfield(struct,{indexList(k)},variable);
        theData = [theData; x];
    end
    if strcmp(variable,'timeInSeconds')
        theData = 1:length(theData);
        theData = theData';
    end
elseif strcmp(operationType,'高搁置')
    indices = find(strcmp({struct.type}, '搁置') == 1);
    theData = getfield(struct,{indices(17 - 2*start)},variable);
elseif strcmp(operationType,'低搁置')
    indices = find(strcmp({struct.type}, '搁置') == 1);
    theData = getfield(struct,{indices(18 - 2*start)},variable);
elseif ~ismember(existType, operationType)
    error(sprintf('The given type does not exist: %s...',operationType));
end



end