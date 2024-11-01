% 
%  input: struct, a struct like 'data', which contains dynamic test data.
%       : start,  2~8, 2 for 20%, 8 for 80%
%       : operationType: '模拟工况' '恒流充电' '恒功率充电' '搁置'
% output: index, 第start个开始的索引
%       : indexList, 第start个开始的连续索引

function [index, indexList] = findPosition(struct,start,operationType)

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
elseif strcmp(operationType,'搁置')
    
elseif ~ismember(existType, operationType)
    error(sprintf('The given type does not exist: %s...',operationType));
end

end

