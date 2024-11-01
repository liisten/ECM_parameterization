% (i) identify different ECMs with IGWO algorithm, on dynamic discharge test.
% These ECM models include,
% 1, average OCV model + 2RC
% 2, discharge OCV model + 2RC
% 4, one state model + 2RC, hysteresis and RC are identified separately
% 5, one state model + 2RC, all params identified at once
% 6, average OCV model + 2RC, try model best
% 7, discharge OCV model + 2RC, try model best
% PS, model no.5 is clearly overfitting, run PSSoverfit.m to see.
% (ii) plot modeling results

clear;
addpath(genpath('E:\Desktop\ECM_parameterization'));

% load the model
temperatureList = [0, 10, 25, 40, 55];
temperature = temperatureList(3); % choose a temperature from the list
cellNum = 1; % choose a cell, 1~4. 1,2: 18650; 3,4: 26650
load(sprintf('./modelData/modelC%d.mat',cellNum)); % load model of the chosen cell

data = load("data.mat");
global model volt curr
volt = data.volt;
curr = data.curr2; % discharge is positive, charge is negative.
soc = data.soc;

Q = 0.98; % cell capacity
timeStep = 1; % the sampling time is 1 s

%% optimize! using multi-step method
% multi-step: first hysteresis related, then RC pairs.
% Number of search agents
N_agts = 50;
% Name of the test function
funcName={'F1','F2','F3','F4'};
% Maximum number of iterations
maxIteration = 500;
% each function corresponding to a different optimization objective
for k = 1:4  
    [lb,ub,dim,fobj] = getECMdetail(funcName{k});
    tic;
    [Fbest,Lbest,Convergence_curve] = IGWOforECM(dim,N_agts,maxIteration,lb,ub,fobj, ...
        temperature,timeStep,volt,curr,model,soc(1),zeros(2,1));
    disp(['The best solution obtained by I-GWO is : ', num2str(Lbest)]);
    disp(['The best optimal value of the objective funciton found by I-GWO is : ',num2str(Fbest)]);
    elapsedTime = toc;
    % save the optimization results
    result_file_path = sprintf('optParams_2/%sT%dN%diter%d', ...
        funcName{k},temperature,N_agts,maxIteration);
    save(result_file_path,'Lbest','Fbest','Convergence_curve', ...
        'elapsedTime');
end
%% simulate ECM with the identified params
% Lbest:(1)R0, (2)R1, (3)C1, (4)R2, (5)C2
load('./optParams_2/F1T25N50iter500.mat');
opt.OCVmodel = 'average';
opt.hysIdtf = 'separate';
[V_term_1, OCV_1] = simECM(Lbest,temperature,timeStep,curr,model,soc(1),zeros(2,1),opt);

opt.OCVmodel = 'discharge';
opt.hysIdtf = 'separate';
[V_term_2, OCV_2] = simECM(Lbest,temperature,timeStep,curr,model,soc(1),zeros(2,1),opt);

load('./optParams_2/F5T25N50iter500.mat');
opt.OCVmodel = 'PSS';
opt.hysIdtf = 'separate';
[V_term_3, OCV_3] = simECM(Lbest,temperature,timeStep,curr,model,soc(1),zeros(2,1),opt);

load('./optParams_2/F3T25N50iter500.mat');
opt.OCVmodel = 'PSS';
opt.hysIdtf = 'together';
[V_term_4, OCV_4] = simECM(Lbest,temperature,timeStep,curr,model,soc(1),zeros(2,1),opt);

load('./optParams_2/F2T25N50iter500.mat');
opt.OCVmodel = 'average';
opt.hysIdtf = 'separate';
[V_term_5, OCV_5] = simECM(Lbest,temperature,timeStep,curr,model,soc(1),zeros(2,1),opt);

load('./optParams_2/F4T25N50iter500.mat');
opt.OCVmodel = 'discharge';
opt.hysIdtf = 'separate';
[V_term_6, OCV_6] = simECM(Lbest,temperature,timeStep,curr,model,soc(1),zeros(2,1),opt);

save("./plot/results_model_1.mat");
%% calculate metrics
V_err_1 = abs(V_term_1-volt) * 1000;
V_err_2 = abs(V_term_2-volt) * 1000;
V_err_3 = abs(V_term_3-volt) * 1000;
V_err_4 = abs(V_term_4-volt) * 1000;
V_err_5 = abs(V_term_5-volt) * 1000;
V_err_6 = abs(V_term_6-volt) * 1000;


% mre 
m_avg_MRE = sum(V_err_1./volt)/length(V_err_1);
m_dis_MRE = sum(V_err_2./volt)/length(V_err_2);
m_pss_MRE = sum(V_err_3./volt)/length(V_err_3);
m_pss_all_MRE = sum(V_err_4./volt)/length(V_err_4);
m_avg_b_MRE = sum(V_err_5./volt)/length(V_err_5);
m_dis_b_MRE = sum(V_err_6./volt)/length(V_err_6);

% rmse in mV
m_avg_RMSE = calc_RMSE(V_term_1,volt) *1000;
m_dis_RMSE = calc_RMSE(V_term_2,volt)*1000;
m_pss_RMSE = calc_RMSE(V_term_3,volt)*1000;
m_pss_all_RMSE = calc_RMSE(V_term_4,volt)*1000;
m_avg_b_RMSE = calc_RMSE(V_term_5,volt)*1000;
m_dis_b_RMSE = calc_RMSE(V_term_6,volt)*1000;

%% save RC parameters to struct model

load('./optParams_2/F1T25N50iter500.mat');
R0 = Lbest(1); R1 = Lbest(2); C1 = Lbest(3);
R2 = Lbest(4); C2 = Lbest(5);
tau1 = exp(-timeStep/R1/C1);
tau2 = exp(-timeStep/R2/C2);

model.(['T', num2str(temperature)]).R0Param = R0;
model.(['T', num2str(temperature)]).tau1Param = tau1;
model.(['T', num2str(temperature)]).tau2Param = tau2;
model.(['T', num2str(temperature)]).R1Param = R1;
model.(['T', num2str(temperature)]).R2Param = R2;
model.(['T', num2str(temperature)]).C1Param = C1;
model.(['T', num2str(temperature)]).C2Param = C2;

modelName = ['modelData/','modelC', num2str(cellNum), '.mat'];
save(modelName, 'model');
fprintf('%s has been saved! \n',modelName);
