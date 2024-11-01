% output a OCV for a given SOC value using the discharge OCV
% input: soc, given soc value, scalar or vector.
%      : temp, scalar, temperature from 0 to 60 degree
%      : model, a struct that contains battery model parameters
% output: corresponding OCV value.

function ocv = OCVfromSOCtemp_DIS(soc,temp,model)
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
ocv = zeros(size(input_soc));
SOC = model.(['T',num2str(temp)]).SOC;
OCV = model.(['T',num2str(temp)]).OCVd; % discharge OCV

I1 = find(input_soc <= SOC(1));
I2 = find(input_soc >= SOC(end));
I3 = find(input_soc > SOC(1) & input_soc < SOC(end));

% for voltages less than lowest stored soc datapoint, extrapolate off 
dsoc = SOC(2) - SOC(1);
if ~isempty(I1)
    dv = OCV(2) - OCV(1);
    ocv(I1) = (input_soc(I1)-SOC(1)).*dv/dsoc + OCV(1);
end

% for voltages greater than highest stored soc datapoint, extrapolate off
if ~isempty(I2)
    dv = OCV(end) - OCV(end-1);
    ocv(I2) = (input_soc(I2)-SOC(end)).*dv/dsoc + OCV(end);
end

% for normal soc range
% if isscalar(input_soc)
%     [~, idx] = sort(abs(SOC - input_soc));
%     idx = idx(1:2); % the first 2
%     idx = sort(idx);
%     p1 = (input_soc - SOC(idx(1)))/dsoc;
%     p2 = (SOC(idx(2)) - input_soc)/dsoc;
%     ocv = disOCV(idx(1))*(1-p1) + disOCV(idx(2))*(1-p2);
% end

if isscalar(input_soc)
    [~, idx] = sort(abs(SOC - input_soc));
    idx = idx(1:2);
    soc1 = SOC(idx(1)); soc2 = SOC(idx(2));
    ocv1 = OCV(idx(1)); ocv2 = OCV(idx(2));

    if soc1 > soc2      
        [soc1, soc2] = deal(soc2, soc1);
        [ocv1, ocv2] = deal(ocv2, ocv1);
    end

    if input_soc > soc1 && input_soc < soc2
        ocv = ocv1 + (input_soc - soc1) * (ocv2 - ocv1) / (soc2 - soc1);
    elseif abs(input_soc - soc1) < 1e-4
        ocv = ocv1;
    elseif abs(input_soc - soc2) < 1e-4
        ocv = ocv2;
    else
        ocv = mean([ocv1,ocv2]);
    end
elseif isvector(input_soc)
    for k=1:length(input_soc)
        [~, idx] = sort(abs(SOC - input_soc(k)));
        idx = idx(1:2); % the first 2
        soc1 = SOC(idx(1)); soc2 = SOC(idx(2));
        ocv1 = OCV(idx(1)); ocv2 = OCV(idx(2));

        if soc1 > soc2
            [soc1, soc2] = deal(soc2, soc1);
            [ocv1, ocv2] = deal(ocv2, ocv1);
        end

        if input_soc(k) > soc1 && input_soc(k) < soc2
            ocv(k) = ocv1 + (input_soc(k) - soc1) * (ocv2 - ocv1) / (soc2 - soc1);
        elseif abs(input_soc(k) - soc1) < 1e-4
            ocv(k) = ocv1;
        elseif abs(input_soc(k) - soc2) < 1e-4
            ocv(k) = ocv2;
        else
            ocv(k) = mean([ocv1,ocv2]);
        end
    end
end

% check if the output ocv value is reasonable
if any(ocv < 1.9 | ocv > 3.7)
    error('输出的OCV范围有误');
end
end

