% find the corresponding SOC for a given OCV value using the average OCV
% where average model is the average of charge and discharge OCV
%  input: ocv, a scalar or a vector 
%       : temp, scalar, from 0 to 60 degree
%       : model, a struct that contains battery model parameters
% output: soc, corresponding SOC value, with identical shape of ocv.

function soc = SOCfromOCVtemp_AVG(ocv,temp,model)
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
% idx = find(model.Temperature == temp,1);
input_ocv = ocv(:);
soc = zeros(size(input_ocv));
SOC = model.(['T',num2str(temp)]).SOC;
OCV = model.(['T',num2str(temp)]).OCVavg;

I1 = find(input_ocv <= OCV(1));
I2 = find(input_ocv >= OCV(end));
I3 = find(input_ocv > OCV(1) & input_ocv < OCV(end));

% for ocvs lower than lowest voltage
docv = OCV(2) - OCV(1);
if ~isempty(I1)
    dz = SOC(2) - SOC(1);
    soc(I1) = (input_ocv(I1)-OCV(1)).*dz/docv + SOC(1);
end

% for ocvs higher than highest voltage
if ~isempty(I2)
    dz = SOC(end)-SOC(end-1);
    soc(I2) = (input_ocv(I2)-OCV(end)).*dz/docv + SOC(end);
end

% for normal ocv range, manually interpolate
% I4 = (input_ocv(I3) - OCVavg(1)) / docv;
% I5 = floor(I4); I45 = I4 - I5; omI45 = 1 - I45;
% soc(I3) = SOC(I5+1) .* omI45 + SOC(I5+2) .* I45;

% simple but not accurate enough
% if isscalar(input_ocv)
%     [~, idx] = sort(abs(OCVavg - input_ocv));
%     idx = idx(1:2); % the first 2
%     soc = mean(SOC(idx));
% end

% manually interpolate, 
if isscalar(input_ocv)
    [~, idx] = sort(abs(OCV - input_ocv));
    idx = idx(1:2); % the first 2
    ocv1 = OCV(idx(1)); ocv2 = OCV(idx(2));
    soc1 = SOC(idx(1)); soc2 = SOC(idx(2));

    if ocv1 > ocv2
        [ocv1, ocv2] = deal(ocv2, ocv1);
        [soc1, soc2] = deal(soc2, soc1);
    end

    % if abs(ocv1 - ocv2) < 1e-4
    if input_ocv > ocv1 && input_ocv < ocv2
        soc = soc1 + (input_ocv - ocv1) * (soc2 - soc1) / (ocv2 - ocv1);
    else
        soc = mean([soc1,soc2]);    
    end
elseif isvector(input_ocv)
    for k = 1:length(input_ocv)
        [~, idx] = sort(abs(OCV - input_ocv(k)));
        idx = idx(1:2); % the first 2
        ocv1 = OCV(idx(1)); ocv2 = OCV(idx(2));
        soc1 = SOC(idx(1)); soc2 = SOC(idx(2));

        if ocv1 > ocv2
            [ocv1, ocv2] = deal(ocv2, ocv1);
            [soc1, soc2] = deal(soc2, soc1);
        end

        if input_ocv(k) > ocv1 && input_ocv(k) < ocv2
            soc(k) = soc1 + (input_ocv(k) - ocv1) * (soc2 - soc1) / (ocv2 - ocv1);
        else
            soc(k) = mean([soc1,soc2]);
        end
    end
end

% check if the output soc value is reasonable
if any(soc < -0.1 | soc > 1.1)
    error('输出的SOC范围有误');
end

end