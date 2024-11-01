% This function containts full information and implementations of the benchmark
% lb is the lower bound: lb=[lb_1,lb_2,...,lb_d]
% up is the uppper bound: ub=[ub_1,ub_2,...,ub_d]
% dim is the number of variables (dimension of the problem)

function [lb,ub,dim,fobj] = getECMdetail(funcName)

switch funcName
    case 'F1'
        fobj = @F1;
        % R0,R1,C1,R2,C2
        lb = [0,0,100,0,100];
        ub = [0.05,0.05,5000,0.05,5000];
        dim = 5;
    case 'F2'
        fobj = @F2;
        % R0,R1,C1,R2,C2,G,M,M0
        lb = [0,0,100,0,100,0,0,0];
        ub = [0.05,0.05,5000,0.05,5000,500,0.05,0.05];
        dim = 8;
    case 'F3'
        fobj = @F3;
        % R0,R1,C1,R2,C2
        lb = [0,0,100,0,100];
        ub = [0.05,0.05,5000,0.05,5000];
        dim = 5;
    case 'F4'
        fobj = @F4;
        % R0,R1,C1,R2,C2, with preoptimized G,M,M0
        lb = [0,0,100,0,100];
        ub = [0.05,0.05,5000,0.05,5000];
        dim = 5;
end

end

% F1, 2RC with average OCV
function z = F1(x,temp,deltaT,vk,ik,model,z0,iR0)
R0 = x(1); R1 = x(2); C1 = x(3);
R2 = x(4); C2 = x(5);
ik = ik(:); iR0 = iR0(:);
Q = 0.98;

RCfact = [exp(-deltaT/R1/C1);exp(-deltaT/R2/C2)];
% calculate irk
irk=zeros([length(ik) length(iR0)]); irk(1,:) = iR0;
for k = 2:length(ik)
    irk(k,:) = RCfact'.*irk(k-1,:) + (1-RCfact')*ik(k-1);
end

% calculate soc
soc = z0-cumsum([0;ik(1:end-1)])*deltaT/(Q*3600);
if any(soc>1.1)
    warning('Current may have wrong sign as SOC > 110%');
end

for k=1:length(soc)
    OCV(k,1) = OCVfromSOCtemp_AVG(soc(k),temp,model);
end
vm = OCV - ik.*R0 - irk*[R1;R2];
z = calc_RMSE(vk,vm);
end

% F2, one state model + 2RC, optimize all params at once
function z = F2(x,temp,deltaT,vk,ik,model,z0,iR0)
R0 = x(1); R1 = x(2); C1 = x(3);
R2 = x(4); C2 = x(5);
G = x(6); M = x(7); M0 = x(8);
ik = ik(:); iR0 = iR0(:);
Q = 0.98;

RCfact = [exp(-deltaT/R1/C1);exp(-deltaT/R2/C2)];
% calculate irk
irk=zeros([length(ik) length(iR0)]); irk(1,:) = iR0;
for k = 2:length(ik)
    irk(k,:) = RCfact'.*irk(k-1,:) + (1-RCfact')*ik(k-1);
end

% calculate soc
soc = z0-cumsum([0;ik(1:end-1)])*deltaT/(Q*3600);
if any(soc>1.1)
    warning('Current may have wrong sign as SOC > 110%');
end

H0 = 0;
Hys = zeros([length(ik) 1]); 
Hys(1) = H0;
sI = 0*Hys;
fac = exp(-abs(G*ik*deltaT/(3600*Q)));
for k = 2:length(ik)
    Hys(k) = fac(k-1)*Hys(k-1) - (1-fac(k-1))*sign(ik(k-1));
    sI(k) = sign(ik(k));
    if abs(ik(k))<Q/100
        sI(k) = sI(k-1);
    end
end

for k = 1:length(soc)
    OCV(k,1) = OCVfromSOCtemp_AVG(soc(k),temp,model);
end
vm = OCV - ik * R0 - irk*[R1;R2] + M * Hys + M0 * sI;
z = calc_RMSE(vk,vm);
end

% F3, discharge OCV model + 2RC
function z = F3(x,temp,deltaT,vk,ik,model,z0,iR0)
R0 = x(1); R1 = x(2); C1 = x(3);
R2 = x(4); C2 = x(5);
ik = ik(:); iR0 = iR0(:);
Q = 0.98;

RCfact = [exp(-deltaT/R1/C1);exp(-deltaT/R2/C2)];
% calculate irk
irk=zeros([length(ik) length(iR0)]); irk(1,:) = iR0;
for k = 2:length(ik)
    irk(k,:) = RCfact'.*irk(k-1,:) + (1-RCfact')*ik(k-1);
end

% calculate soc
soc = z0-cumsum([0;ik(1:end-1)])*deltaT/(Q*3600);
if any(soc>1.1)
    warning('Current may have wrong sign as SOC > 110%');
end

for k=1:length(soc)
    OCV(k,1) = OCVfromSOCtemp_DIS(soc(k),temp,model);
end
vm = OCV - ik.*R0 - irk*[R1;R2];
z = calc_RMSE(vk,vm);
end

% F4, one state model + 2RC
% with preoptimized one state model parameters
% and only optimize RC params here
function z = F4(x,temp,deltaT,vk,ik,model,z0,iR0)
R0 = x(1); R1 = x(2); C1 = x(3);
R2 = x(4); C2 = x(5);
G = model.(['T',num2str(temp)]).GParam; 
M = model.(['T',num2str(temp)]).MParam; 
M0 = model.(['T',num2str(temp)]).M0Param;
ik = ik(:); iR0 = iR0(:);
Q = 0.98;

RCfact = [exp(-deltaT/R1/C1);exp(-deltaT/R2/C2)];
% calculate irk
irk=zeros([length(ik) length(iR0)]); irk(1,:) = iR0;
for k = 2:length(ik)
    irk(k,:) = RCfact'.*irk(k-1,:) + (1-RCfact')*ik(k-1);
end

% calculate soc
soc = z0-cumsum([0;ik(1:end-1)])*deltaT/(Q*3600);
if any(soc>1.1)
    warning('Current may have wrong sign as SOC > 110%');
end

H0 = 0;
Hys = zeros([length(ik) 1]); 
Hys(1) = H0;
sI = 0*Hys;
fac = exp(-abs(G*ik*deltaT/(3600*Q)));
for k = 2:length(ik)
    Hys(k) = fac(k-1)*Hys(k-1) - (1-fac(k-1))*sign(ik(k-1));
    sI(k) = sign(ik(k));
    if abs(ik(k))<Q/100
        sI(k) = sI(k-1);
    end
end

for k = 1:length(soc)
    OCV(k,1) = OCVfromSOCtemp_AVG(soc(k),temp,model);
end
vm = OCV - ik * R0 - irk*[R1;R2] + M * Hys + M0 * sI;
z = calc_RMSE(vk,vm);
end