clear all; clc;

gamma = 0.9;
beta = 0.49;
xi = 0.3;

dt_T = 0.1;
Omega = dt_T * 2 * pi;

[eigen_1, eigen_2] = newmark_eigen(beta, gamma, Omega, xi)

[ A ] = newmark_A(beta, gamma, Omega, xi);


E = eig(A)
