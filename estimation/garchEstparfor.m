function [thetahat, hhat, sumLLF, epsihat] = garchEstparfor(data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: 
%   This function estimate a GARCH(1,1)  process under constraints. It is
%   coded in order to use parallel computing in MATLAB.
%   If needed, the user can use the robustvcv.m functon of the MFET toolbox
%   coded by Kevin Sheppard (link to the toolbox: )
%
% INPUTS:
%   data: The data of the process
%
% OUTPUTS:
%   thetahat: Vector of estimate coefficient
%   hhat: Vector composed of estimated condition volatility 
%   sumLLF: Scalair of the log-likelihood value
%   epsihat: Vector of estimated data
%
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-fcomte.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check the inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global y ;
global ORDERS;
global nbcoefgarch;

y = data;
ORDERS = [1 1];
nbcoefgarch = sum(ORDERS);
valdep = [0.01 0.1 0.7]';

LB  = zeros(3,1);
UB  = ones(3,1);
% We call fmincon function for the estimation
options  =  optimset('fmincon');
options  =  optimset(options , 'Algorithm ','interior-point');
options  =  optimset(options , 'TolFun'      , 1e-006);
options  =  optimset(options , 'TolX'        , 1e-006);
options  =  optimset(options , 'TolCon'      , 1e-006);
options  =  optimset(options , 'Display'     , 'off');
options  =  optimset(options , 'Diagnostics' , 'off');
options  =  optimset(options , 'LargeScale'  , 'off');
options  =  optimset(options , 'MaxIter'     , 1500);
options  =  optimset(options , 'Jacobian'     ,'off');
options  =  optimset(options , 'MeritFunction'     ,'multiobj');
options  =  optimset(options , 'MaxFunEvals' , 3000);
A = [-eye(3) ; zeros(1,1) ones(1,2)];
bt = [zeros(1, 3) 1 - (1e-6)];
LB  = zeros(3,1);
UB  = [];
thetahat = fmincon('garchLik',valdep,A,bt,[],[],LB,UB,[],options); 
[sumLLF,~,hhat,epsihat] = garchLik(thetahat);
%[~,~,~,scores,~,gross_scores]=robustvcv('garchLik',thetahat,0);

end

