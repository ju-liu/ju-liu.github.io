clear all; clc;
syms af am dt gamma lambda

A = [ -lambda * af, am / dt; 1, -gamma ];

B = [ lambda * (1-af), (am-1)/dt; 1, 1-gamma];


d = det(A);


C = simplify( inv(A) * B );


trace(C) * d
