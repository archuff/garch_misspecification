function [simulatedata, ht] = garch_simulate_rsample(t,parameters,etah,p,q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: 
%   This function simulates GARCH(p,q) model with resampling from estimated residuals.
%   It is an auxilary function used for bootstrap procedure. 
%
%
% INPUTS:
%   t: length of the process simulated
%   parameters: Vecotr of the paramaters 
%   etah: Vector of residuals estimated
%   p, q: orders or GARCH process (scalars)
%
% OUTPUTS:
%   simulatedata: Vector of simulated data with GARCH variance
%   ht: Vector of the conditional variance 
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-fcomte.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Separate the parameters
omega=parameters(1);
alpha=parameters(2:p+1);
beta=parameters(p+2:p+q+1);

%Initialize the random numbers
m  =  max([p,q]);
RandomNums = datasample(etah,t+m);
UncondStd =  sqrt(omega/(1-sum(alpha)-sum(beta)));
h=UncondStd.^2*ones(t+m,1);
e=UncondStd*ones(t+m,1);
T=size(e,1);
parameters=[omega;alpha;beta];
for t = (m + 1):T
    h(t) = parameters' * [1 ; e(t-(1:p)).^2 ; h(t-(1:q))];
    e(t)= RandomNums(t)*sqrt(h(t));
end
simulatedata=e((m+1):T);
ht=h(m+1:T);