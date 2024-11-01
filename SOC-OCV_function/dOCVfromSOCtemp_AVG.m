% compute an approximation to the derivative of the OCV to SOC
% inputs: soc, given soc value, scalar or vector
%       : temp, scalar, temperature from 0 to 60 degree
%       : model, a struct that contains battery model parameters
% output: corresponding dOCV value.

function docv = dOCVfromSOCtemp_AVG(soc,temp,model)
if ~isscalar(temp)
    error('temp has to be a scalar!');
end
if temp >=0 && temp < 10
    temp = 0;
elseif temp >= 10 && temp < 25
    temp = 10;
elseif temp >= 25 && temp < 40
    temp = 25;
elseif temp >= 40 && temp < 55
    temp = 40;
elseif temp >= 55 && temp < 60
    temp = 55;
else
    error('temp has to be within 0 and 60!')
end

input_soc = soc(:);
docv = zeros(size(input_soc));
SOC = model.(['T',num2str(temp)]).SOC;
dOCV = model.(['T',num2str(temp)]).dOCV;

I1 = find(input_soc <= SOC(1));
I2 = find(input_soc >= SOC(end));
I3 = find(input_soc > SOC(1) & input_soc < SOC(end));

% for voltages less than lowest stored soc datapoint, extrapolate off 
dsoc = SOC(2) - SOC(1);
if ~isempty(I1)
    docv(I1) = docv(1);
end

% for voltages greater than highest stored soc datapoint, extrapolate off
if ~isempty(I2)
    docv = docv(end);
end

% for normal soc range
% if isscalar(input_soc)
%     [~, idx] = sort(abs(SOC - input_soc));
%     idx = idx(1:2); % the first 2
%     idx = sort(idx);
%     p1 = (input_soc - SOC(idx(1)))/dsoc;
%     p2 = (SOC(idx(2)) - input_soc)/dsoc;
%     docv = dOCV(idx(1))*(1-p1) + dOCV(idx(2))*(1-p2);
% end

if isscalar(input_soc)
    [~, idx] = sort(abs(SOC - input_soc));
    idx = idx(1:2);
    soc1 = SOC(idx(1)); soc2 = SOC(idx(2));
    docv1 = dOCV(idx(1)); docv2 = dOCV(idx(2));

    if soc1 > soc2      
        [soc1, soc2] = deal(soc2, soc1);
        [docv1, docv2] = deal(docv2, docv1);
    end

    if input_soc > soc1 && input_soc < soc2
        docv = docv1 + (input_soc - soc1) * (docv2 - docv1) / (soc2 - soc1);
    elseif abs(input_soc - soc1) < 1e-4
        docv = docv1;
    elseif abs(input_soc - soc2) < 1e-4
        docv = docv2;
    else
        docv = mean([docv1,docv2]);
    end
elseif isvector(input_soc)
    for k=1:length(input_soc)
        [~, idx] = sort(abs(SOC - input_soc(k)));
        idx = idx(1:2); % the first 2
        soc1 = SOC(idx(1)); soc2 = SOC(idx(2));
        docv1 = dOCV(idx(1)); docv2 = dOCV(idx(2));

        if soc1 > soc2
            [soc1, soc2] = deal(soc2, soc1);
            [docv1, docv2] = deal(docv2, docv1);
        end

        if input_soc(k) > soc1 && input_soc(k) < soc2
            docv(k) = docv1 + (input_soc(k) - soc1) * (docv2 - docv1) / (soc2 - soc1);
        elseif abs(input_soc(k) - soc1) < 1e-4
            docv(k) = docv1;
        elseif abs(input_soc(k) - soc2) < 1e-4
            docv(k) = docv2;
        else
            docv(k) = mean([docv1,docv2]);
        end
    end



end

