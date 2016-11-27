function [ U ] = encodeU( rho, A, V, T , dA)
%ENCODE Summary of this function goes here
%   Detailed explanation goes here

gamma=1.4;

U1=rho.*A;

U2=rho.*A.*V;

U3=rho.*(T./(gamma-1)+gamma/2*V.^2).*A;

U=[U1;U2;U3];


end

