function [sumLLF, LLF, h, epsihat] = garchLik(para)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: 
%   This function computes the likelihood of a GARCH(1,1) process. This
%   functin is coded in order to be used in parallel computing in MATLAB.
%
% INPUTS:
%   para: vector of the parameters
%
% OUTPUTS:
%   sumLLF: The sum of the individual likelihood
%   LLF: The likelihood in theta 
%   h: h estimated
%   epsihat: esimated residuals
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-fcomte.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global y
global ORDERS
global nbcoefgarch

dim = length(y);
maxi = max(ORDERS);

h = zeros(length(y),1);
h(1:maxi) = var(y);
mu = zeros(dim,1);
p = ORDERS(1);
q = ORDERS(2);
for t = (maxi + 1):dim,
	h(t) = para' * [1 ; (y(t-(1:p))).^2 ; h(t-(1:q))];
end

t = 1:dim;
LLF = 0.5*((log(h(t))) + ((y(t).^2)./h(t)) + log(2*pi));
sumLLF = (sum(LLF)); 
if isnan(LLF)
    sumLLF=1e+6;
end

epsihat = (y(t))'; 

end