function [ x ] = combmj(m,j)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes combination between two real (j among m). 
%
%
%   INPUT:
%       - m the number of element in the set
%       - j: the number of distinct element in the S
%
%   OUTPUT:
%       - x, the x-combination possible
%
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = factorial(m)/(factorial(j)*factorial(m-j));

end

