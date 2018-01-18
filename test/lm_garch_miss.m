function results = lm_garch_miss(e,h,param,ORDTEST,ORDTAY,k,flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: 
%   This function computes the statitiscal test of misspecification 
%   in GARCH-type model by Chuffart, Flachaire and Peguin-Feissolle.
%       H0: The error term follows a GARCH(p,q) process
%       H1: The conditional variance is misspecified
%   This function uses two auxiliary functions: chuflapeg_core and
%   taylor2var and a function from MATLAB: pca
%
% INPUTS:
%   e: Vector of estimated residuals
%   h: Vector of estimated condition volatility under the null hypothesis
%   param: Vector of estimated parameters under the null hypothesis
%   ORDTEST: Vector of the lags of the null hypothesis
%   ORDTAY: Vector of the lags in the Taylor expansion 
%   k: Sclalar of the Taylor expansion order
%   flag: Boolean (1 or 0), indicated if the PCA has to be used
%
% OUTPUTS:
%   results: a structure with different components like the value of the
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-fcomte.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ptest = ORDTEST(2);
qtest = ORDTEST(1);
ptay = ORDTAY(2);
qtay = ORDTAY(1);
nblag = ptest+qtest+1;
e2 = e.^2;
T = length(e2);
elagtay = zeros(T,qtay);
hlagtay = zeros(T,ptay);
elagtest = zeros(T,qtest);
hlagtest = zeros(T,ptest);

for t = qtay+1:T,
    elagtest(t,:) = e(t-(1:qtest))';
    hlagtest(t,:) = h(t-(1:ptest))';
end
for t = qtay+1:T,
    elagtay(t,:) = e(t-(1:qtay))';
    hlagtay(t,:) = h(t-(1:ptay))';  
end

e2lagtest = elagtest.^2;
to = [e2lagtest hlagtest];
cto = size(to,2);
xx = [elagtay hlagtay];
taylormatrix = taylor2var(xx,xx,k);
v1 = taylormatrix;
for i = 1:cto,
    cv = size(v1,2);
    for j = 1:cv
        if to(:,i) == v1(:,j)
            v1(:,j) = [];
            break
        end
    end
end

if flag == 1,
   [~,score,latent] = pca(v1);
   r=find(cumsum(latent)./sum(latent) > 0.90); 
   prinCom = score(:,1:r);
   dof = size(prinCom,2);
else
    prinCom = v1;
    dof = size(prinCom,2);
end

sres = e2./h - 1;
dgarch = [ones(T,1) e2 h];
stat = lm_garch_miss_core(T,h,param,prinCom,dgarch,sres,nblag) ;
pval = 1-chi2cdf(stat,dof);
results.pval = pval;
results.stat = stat;
results.v1 = v1;
results.d1 = results.pval < 0.01;
results.d5 = results.pval < 0.05;
results.d10 = results.pval < 0.1;   
results.prinCom = prinCom;
results.dof = dof;


end
