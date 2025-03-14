% (1) identify RC parameters of different ECMs with IGWO algorithm
% only identify RC on a single condition (sc). 
% i.e. one temperature, one usage condition
% using dynamic discharge profiles for identification
% These ECM models include,
% 1, average OCV model + 2RC
% 2, discharge OCV model + 2RC
% 3, API OCV model + 2RC (the proposed)
% 4, CPI OCV model + 2RC
% 5, GPI OCV model + 2RC
% 6, one state model + 2RC, hysteresis and RC are identified separately

clear;
addpath(genpath('../Code'));

% load the model
temp_list = [0, 10, 25, 40, 55];
cell_list = {'cell-1-5','cell-1-10','cell-2-5','cell-2-10'};
cell_num = 1; % choose a cell, 1~4. 1,2: 18650; 3,4: 26650
load(sprintf('./modelData/modelC%d_pi_1.mat',cell_num)); % load model with PI OCV submodel

% conditions
main_dir_1 = 'D:\迟滞实验\迟滞+SOC估计测试\变充电+US06放电\处理后数据';
main_dir_2 = 'D:\迟滞实验\迟滞+SOC估计测试\恒流充电+变放电\处理后数据';
cond_chg = {'CC-1C','CC-2C','CC-3C','CP-3or6W','CP-6or12W','CP-9or18W','MCC','MCP'};
cond_dis = {'CLTC','NEDC','US06','WLTC'};
temp = temp_list(3); % choose a temperature
%% run this if optimize at a specific condtion
% load data

cond_chg_num = 1; % choose a scenario, 1~8
chg_file_path = fullfile(main_dir_1,num2str(temp), ...
    cond_chg{cond_chg_num},[cell_list{cell_num},'.mat']);
if ~exist(chg_file_path,"file")
    disp(['mat file does not exist: ', chg_file_path]);
end
load(chg_file_path); % load the struct 'data' of dynamic test

% discharge data to 20%
raw_curr = getDynData(data,'current',2,'模拟工况');
curr = - raw_curr; % make discharge current to be positive (+)
volt = getDynData(data,'voltage',2,'模拟工况');
soc = getDynData(data,'soc',2,'模拟工况');
Q = model.(['T',num2str(temp)]).QParam; % cell capacity
time_step = 1; % the sampling time is 1 s

% optimize! using multi-step method
% multi-step: first hysteresis related, then RC pairs.

N_agts = 50; % Number of search agents
funcName={'F1','F2','F3','F4','F5','F6'}; % Name of the test function
maxIteration = 500; % Maximum number of iterations

% each function corresponding to a different optimization objective
for k = 1:6
    [lb,ub,dim,fobj] = getECMdetail(funcName{k});
    tic;
    [Fbest,Lbest,Convergence_curve] = IGWOforECM(dim,N_agts,maxIteration,lb,ub,fobj, ...
        temp,time_step,volt,curr,model,soc(1),zeros(2,1));
    disp(['The best solution obtained by I-GWO is : ', num2str(Lbest)]);
    disp(['The best optimal value of the objective funciton found by I-GWO is : ',num2str(Fbest)]);
    elapsed_time = toc;
    % save the optimization results
    result_file_path = sprintf('optParams/%sT%dN%diter%d', ...
        funcName{k},temp,N_agts,maxIteration);
    save(result_file_path,'Lbest','Fbest','Convergence_curve', ...
        'elapsed_time');
end
%% save RC parameters to struct model

load('./optParams/F1T25N50iter500.mat');
R0 = Lbest(1); R1 = Lbest(2); C1 = Lbest(3);
R2 = Lbest(4); C2 = Lbest(5);
tau1 = exp(-time_step/R1/C1);
tau2 = exp(-time_step/R2/C2);

model.(['T', num2str(temp)]).R0Param = R0;
model.(['T', num2str(temp)]).tau1Param = tau1;
model.(['T', num2str(temp)]).tau2Param = tau2;
model.(['T', num2str(temp)]).R1Param = R1;
model.(['T', num2str(temp)]).R2Param = R2;
model.(['T', num2str(temp)]).C1Param = C1;
model.(['T', num2str(temp)]).C2Param = C2;

modelName = ['modelData/','modelC', num2str(cell_num), '_pi_1.mat'];
save(modelName, 'model');
fprintf('%s has been saved! \n',modelName);
