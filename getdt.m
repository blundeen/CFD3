function [ dt ] = getdt( T, V, dx, CFL )
%GETDT Summary of this function goes here
%   Detailed explanation goes here

dt=min(abs(CFL*dx./(sqrt(T)+V)));




end

