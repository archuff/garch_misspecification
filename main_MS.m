% This script replicates the power results in Chuffart et al. when the true
% DGP is a Markov-Switching GARCH model (figure 2). More precisely, it
% replicates one point of the Figure 2.
% To replicate all the points, the user has to change the values of the transition matrix
% P in order to increase or decrease the probability to be in regime 1. 


mycluster = parcluster('local');
delete(mycluster.jobs);
parpool(mycluster, 'open');

clear all

T = 1000;
repli = 1000;
nboot = 499;

%Initialization 
NameProcess = 'MS_05m_2';
LM112 = zeros(1,repli);
LM113 = zeros(1,repli);
LM223 = zeros(1,repli);
lt = zeros(1,repli);
ho = zeros(1,repli);
LM112b = zeros(1,nboot);
LM113b = zeros(1,nboot);
LM223b = zeros(1,nboot);
ltb = zeros(1,nboot);
hob = zeros(1,nboot);
asl = zeros(5,nboot);
k = 2;
flag = 2;
params = [0.01 0.1 0.7 ; 0.0001 0.2 0.4];
ORDERS = [0 1 1];
M = [0.95  0.05; 0.05 0.95];
ASL = zeros(5,repli);
for i = 1:repli,
    i 
    [data, hhat, statesHaas, Yhaas] = msGarchSim(T,params,ORDERS,M,k,flag,0);    
    [thetahat, hh, sumLLF, epsihat] = garchEstparfor(data);
    etahat =(data./sqrt(hh));
    LM112Model(i) = lm_garch_miss(data,hh,thetahat,[1 1],[1 1],2,0);            
    LM113Model(i) = lm_garch_miss(data,hh,thetahat,[1 1],[1 1],3,1);
    LM223Model(i) = lm_garch_miss(data,hh,thetahat,[1 1],[2 2],3,1);    
    lt_s(i) = lt_test(data,hh,thetahat,[1 1]); 
    ho_s(i) = ho_test(data,hh,[1 1],thetahat);
    lt(i) = lt_s(i).stat;
    ho(i) = ho_s(i).stat;
    LM112(i) = LM112Model(i).stat;
    LM113(i) = LM113Model(i).stat;
    LM223(i) = LM223Model(i).stat; 
    parfor j = 1:nboot, 
        [dataBoot, varBoot] = garch_simulate_rsample(T, thetahat,etahat, 1, 1);
        [thetahatBoothat, varBoothat, sumLLF, dataBoothat] = garchEstparfor(dataBoot);
        LM112ModelBoot = lm_garch_miss(dataBoot,varBoothat,thetahatBoothat,[1 1],[1 1],2,0); 
        LM113ModelBoot = lm_garch_miss(dataBoot,varBoothat,thetahatBoothat,[1 1],[1 1],3,1); 
        LM223ModelBoot = lm_garch_miss(dataBoot,varBoothat,thetahatBoothat,[1 1],[2 2],3,1); 
        lt_Boot = lt_test(dataBoot,varBoothat,thetahatBoothat,[1 1]);
        ho_Boot = ho_test(dataBoot,varBoothat,[1 1],thetahatBoothat);         
        LM112b(j) = LM112ModelBoot.stat;
        LM113b(j) = LM113ModelBoot.stat;
        LM223b(j) = LM223ModelBoot.stat;
        ltb(j) = lt_Boot.stat;
        hob(j) = ho_Boot.stat;
    end
    asl(1,:) = ltb > lt(i);
    asl(2,:) = hob > ho(i);
    asl(3,:) = LM112b > LM112(i);
    asl(4,:) = LM113b > LM113(i);  
    asl(5,:) = LM223b > LM223(i);    
    ASL(:,i) = sum(asl,2)/nboot;
end

sum(ASL<0.05,2)/repli
save(NameProcess);
