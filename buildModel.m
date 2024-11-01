% import raw data and extract model parameter
% output: "model.mat", contains an struct named "model"
% this struct contains OCV data and identified parameters.
cd('D:\2023_Paper3\Code');
cellNum = 1; % choose a cell, 1~4. 1,2: 18650; 3,4: 26650
temperatureList = [0, 10, 25, 40, 55];
lengthTemp = length(temperatureList);

%% compute capacity and coulombic efficiency first
model = struct();
for k = 1:lengthTemp
    if cellNum == '1'
        model.(['T',num2str(temperatureList(k))]).QParam = 0.98;
    elseif cellNum == '2'
        model.(['T',num2str(temperatureList(k))]).QParam = 0.98;
    elseif cellNum == '3'
        model.(['T',num2str(temperatureList(k))]).QParam = 2.35;
    elseif cellNum == '4'
        model.(['T',num2str(temperatureList(k))]).QParam = 2.35;
    end
    % assume the coulombic efficiency to be 1
    % ce = discharge capacity / charge capacity
    model.(['T',num2str(temperatureList(k))]).etaParam = 1;
end

%% extract charge, discharge, and average OCV
path1 = 'D:\2023_Paper4\FigureProduce\majorLoop';

for k = 1:lengthTemp
    filePath1 = fullfile(path1,sprintf("major-%d-%d.csv", ...
        cellNum, temperatureList(k)));
    % disp(filePath);
    data = readmatrix(filePath1);
    socChgMain = data(:,1);
    ocvChgMain = data(:,2);
    socDisMain = data(:,3);
    ocvDisMain = data(:,4);
    % 确保socChgMain中没有重复元素
    [~, idx] = unique(socChgMain);
    socChgMain = socChgMain(idx);
    ocvChgMain = ocvChgMain(idx);
    socDisMain = socDisMain(idx);
    ocvDisMain = ocvDisMain(idx);

    newSoc = linspace(0,1,1000)';
    newOcvChg = interp1(socChgMain,ocvChgMain,newSoc,'makima');
    newOcvDis = interp1(socDisMain,ocvDisMain,newSoc,'makima');
    params.SOC = newSoc;
    params.OCVc = newOcvChg;
    params.OCVd = newOcvDis;
    params.OCVavg = (newOcvChg + newOcvDis) / 2;
    % calcualte the slope of SOC-OCV curve
    delta_SOC = diff(params.SOC);
    delta_OCV = diff(params.OCVavg);
    docv = delta_OCV./delta_SOC;
    docv = [docv;docv(end)]; % keep the same length with SOC
    UL = 50; % set a upper limit of the maximum slope
    DL = 0.001; % lower limit of the minimum slope
    params.dOCV = max(min(docv,UL),DL);

    model.(['T', num2str(temperatureList(k))]) = params;
    model.Temperature = temperatureList;
end

%% plot three OCV curves 
for k = 1:lengthTemp
    x = model.(['T', num2str(temperatureList(k))]).SOC;
    yc = model.(['T', num2str(temperatureList(k))]).OCVc;
    yd = model.(['T', num2str(temperatureList(k))]).OCVd;
    yavg = model.(['T', num2str(temperatureList(k))]).OCVavg;
    figure(k); hold on; box on;
    plot(x,yc);
    plot(x,yd);
    plot(x,yavg);
    legend('charge','discharge','average',Location='best');
    title(['SOC-OCV curve at ', num2str(temperatureList(k))]);
end
%% load model parameters related to hysteresis: Gamma, M, M0
% please identify these parameters first using hysteresis test data from 
% the same cell! After saving the identified results, we can load them
% from the following folder.
path3 = 'D:\2023_Paper2\Code\PSSParams';

for k = 1:lengthTemp
    mat_files = dir(fullfile(path3,'*.mat'));
    for k1 = 1:length(mat_files)
        file_name = mat_files(k1).name;
        if contains(file_name, sprintf('F1Typechg10T%d', temperatureList(k)))
            full_file_path = fullfile(path3, file_name);
            data = load(full_file_path);
        end
    end
    GParam = data.Lbest(3);
    M0 = data.Lbest(1);
    M = data.Lbest(2);
    % save hysteresis parameters to model
    model.(['T', num2str(temperatureList(k))]).GParam = GParam;
    model.(['T', num2str(temperatureList(k))]).M0Param = M0;
    model.(['T', num2str(temperatureList(k))]).MParam = M;
end

%% model parameters related to RC
% the RC related parameters are identified in idtfRCParam_IGWO.m
% update this structure by running the aforementioned m file.

%% save the model
modelName = ['modelData/','modelC', num2str(cellNum), '.mat'];
save(modelName, 'model');
fprintf('%s has been saved! \n',modelName);