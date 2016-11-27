function [ rho, T, V] = decodeU( U, A )
%DECODEU Summary of this function goes here
%   Detailed explanation goes here
gamma=1.4;

U1=U(1,:);

U2=U(2,:);

U3=U(3,:);

rho=U1./A;
V=U2./(rho.*A);

T=(U3./(rho.*A)-gamma/2*V.^2)*(gamma-1);
end

