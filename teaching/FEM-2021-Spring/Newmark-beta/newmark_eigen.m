function [eigen_1, eigen_2] = newmark_eigen(beta, gamma, Omega, xi)

D = 1 + 2 * gamma * xi * Omega + beta * Omega^2;

A1 = 1 - (xi*Omega + 0.5 * Omega^2 * (gamma + 0.5) ) / D;

A2 = 1 - (2*xi*Omega + Omega^2 * (gamma - 0.5)) / D;

eigen_1 = A1 + sqrt(A1^2 - A2);
eigen_2 = A1 - sqrt(A1^2 - A2);