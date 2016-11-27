function [ F, J ] = encode( U, dA)
%ENCODE Summary of this function goes here
%   Detailed explanation goes here

gamma=1.4;

U1=U(1,:);

U2=U(2,:);

U3=U(3,:);


F1=U2;

F2=U2.^2./U1+(gamma-1)/gamma*(U3-gamma/2*U2.^2./U1);

F3=gamma*U2.*U3./U1-gamma*(gamma-1)/2*U2.^3./U1.^2;

F=[F1;F2;F3];

J1=zeros(1,length(dA));

J2=(gamma-1)/gamma*(U3-gamma/2*U2.^2./U1).*dA;

J3=zeros(1,length(dA));

J=[J1;J2;J3];

end

