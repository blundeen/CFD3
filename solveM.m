function [ res ] = solveM( M )
%SOLVE M Summary of this function goes here
%   Detailed explanation goes here

global A
gamma=1.4;

if M<.00
    res=5;
else
res=A-((gamma+1)/2)^((-gamma-1)/(2*(gamma-1)))*(1+(gamma-1)/2*M^2)^((gamma+1)/(2*(gamma-1)))/M;

end


end

