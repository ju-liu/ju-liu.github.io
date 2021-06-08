clear all; clc;
syms xi Omega beta gamma;
Left  = [ 1, 2*xi * Omega, Omega^2; -beta, 0, 1; -gamma, 1, 0];
Right = [ 0, 0, 0; 0.5-beta, 1,1; 1-gamma, 1, 0];

D = det(Left);

adj_Left = [ -1, Omega^2, 2*xi*Omega; 
    -gamma, gamma*Omega^2, -(1+beta*Omega^2); 
    -beta, -(1+2*gamma*xi*Omega), 2*beta*xi*Omega];

A1 = trace(adj_Left * Right) / D; A1 = simplify(A1);
