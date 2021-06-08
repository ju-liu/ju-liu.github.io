function [A] = newmark_A(beta, gamma, Omega, xi)

Left  = [ 1, 2*xi * Omega, Omega^2; -beta, 0, 1; -gamma, 1, 0];
Right = [ 0, 0, 0; 0.5-beta, 1,1; 1-gamma, 1, 0];

A = inv(Left) * Right;