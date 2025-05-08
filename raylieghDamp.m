function [alpha, beta] = raylieghDamp(omegaN,xi)
% Calculates the Rayleigh damping coeffcient

A = [1, omegaN(1)^2 ; 1 omegaN(2)^2];
b = [2*xi(1)*omegaN(1) ; 2*xi(2)*omegaN(2)];
x = A \ b;
alpha = x(1);
beta = x(2);