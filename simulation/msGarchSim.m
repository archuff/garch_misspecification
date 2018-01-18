function [simulatedEpsi, H, states, Y] = msGarchSim(dim,para,ORDERS,M,k,flag,flagstd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: 
%   This function simulates Markov-Switching Klaassen (2002) 
%   and Haas & al (2004) with Gaussian innovations. This function generates 1000 observations
%   in more for reducing any bias. 
%
% INPUTS:
%   dim: length of the process simulated
%   para: The matrix of the paramaters. Size of Para has to be equal to
%   k*n with k the number of regime, and n the number of GARCH parameter
%   (p+q+1)
%   ORDERS: A vectore of length 3 with the order of the GARCH in the
%   regimes
%   M: The matrice of probability transition of the Markov Chain, has to be
%   of size k*k
%   k: The number of regime
%   flag: To decide what sort type of markov switching Garch model to
%   simulate, flag = 1 for Klaassen and flag = 2 for
%   Haas.
%   s: seed
%
% OUTPUTS:
%   simulateEpsi: Process with MS-GARCH variances
%   H: Conditional variances, a dim*k matrix
%   states: The state of the process at any time a dim*k matrix if haas
%   Y: simulated epis for the moment. ARMA further. 
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-fcomte.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Checkin' INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    c = ORDERS(1);
    p = ORDERS(2);
    q = ORDERS(3);
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    maxi = max(p,q);
    T = dim+10000; 
    states = zeros(T,k); 
    states(1,1) = 1; 
    if flagstd == 1,
        n = stdtrnd(5,T,1);
    else
        n = randn(T,1);
    end
    
    h = zeros(T,k);
    h(1,:) = 1; 
    H = zeros(T,1);
    Y = zeros(T,k);
    simulatedEpsi = zeros(T,1); 

    if c == 1, 
        cons(1:k) = para(1:k,1);
        paraGarch = para(1:k,2:end);
    else
        cons = zeros(1,k);
        paraGarch = para(1:k,1:end);
    end
    switch flag,
        case 1
            states = markovSim(T,M);
            for t = maxi+1:T,
                i = find(states(t,:));
                h(t) = paraGarch(i,:) * [1;simulatedEpsi(t-(1:q)).^2 ;h(t-(1:p))];
                simulatedEpsi(t) = n(t)*sqrt(h(t));
                Y(t) = cons(i) + simulatedEpsi(t);
            end
            H = h(maxi+10000:T,1);
            simulatedEpsi =  simulatedEpsi(maxi+10000:T);
            states = states(maxi+10000:T,:); 
            Y = Y(maxi+10000:T)';
            
        case 2           
            states = markovSim(T,M);        
            for t = maxi+1:T,
                for i = 1:k,
                    h(t,i) = paraGarch(i,:) * [1;simulatedEpsi(t-(1:q))'.^2;h(t-(1:p),i)] ; 
                end
                i = find(states(t,:));
                simulatedEpsi(t) = n(t)*sqrt(h(t,i));
                Y(t) = cons(i) +  simulatedEpsi(t);   
            end
            H =states.*h; 
            for j = 1:T,
                if H(j,1) == 0,
                    H(j,1) = H(j,2);
                end
            end
            H = H(maxi+10000:T,1);
            simulatedEpsi =  simulatedEpsi(maxi+10000:T);
            states = states(maxi+10000:T,:); 
            Y = Y(maxi+10000:T)';
    end
end

