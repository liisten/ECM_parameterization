% simulate ECMs with given parameters 
%
function [vm,OCV] = simECM(paraVec,temp,deltaT,ik,model,z0,iR0,opt)

R0 = paraVec(1); R1 = paraVec(2); C1 = paraVec(3);
R2 = paraVec(4); C2 = paraVec(5);

% the way to obtain hysteresis parameters.
switch opt.hysIdtf
    case 'separate'
        G  = getParaECM('GParam',temp,model);
        M  = getParaECM('MParam',temp,model);
        M0 = getParaECM('M0Param',temp,model);
    case 'together'
        G = paraVec(6); 
        M = paraVec(7); 
        M0 = paraVec(8);
    otherwise
        error('wrong')
end

ik = ik(:); iR0 = iR0(:);

Q = 0.98;

tau1 = exp(-deltaT/R1/C1);
tau2 = exp(-deltaT/R2/C2);
RCfact = [tau1;tau2];

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

% calculate the hysteresis state
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

switch opt.OCVmodel
    case 'average'
        for k= 1:length(soc)
            OCV(k,1) = OCVfromSOCtemp_AVG(soc(k),temp,model);
        end
        vm = OCV - ik * R0 - irk*[R1;R2];
    case 'discharge'
        for k= 1:length(soc)
            OCV(k,1) = OCVfromSOCtemp_DIS(soc(k),temp,model);
        end
        vm = OCV - ik * R0 - irk*[R1;R2];
    case 'PSS'
        for k = 1:length(soc)
            OCV(k,1) = OCVfromSOCtemp_AVG(soc(k),temp,model);
        end
        vm = OCV - ik * R0 - irk*[R1;R2] + M * Hys + M0 * sI;
    otherwise
        error('please choose an OCV model.')
end

end

