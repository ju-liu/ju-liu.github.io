clear all; clc; clf;
hold off;

gamma = 0.5;
beta = [ 0.45, 0.47, 0.49, 0.55];
xi = 0.0;

dt_T = 0.001 : 0.0001 : 0.01;
dt_T = [ dt_T , 0.01 : 0.001 : 0.1];
dt_T = [ dt_T, 0.1 : 0.01 : 1 ];
dt_T = [ dt_T, 1 : 0.1 : 10 ];
dt_T = [ dt_T, 10 : 1 : 100 ];
dt_T = [ dt_T, 100 : 10 : 1000 ];

rho = zeros(length(dt_T),1);

for jj = 1 : length(beta)
    for ii = 1 : length(dt_T)
        Omega = dt_T(ii) * 2 * pi;
        
        [eigen_1, eigen_2] = newmark_eigen(beta(jj), gamma, Omega, xi);
        
        rho(ii) = abs(eigen_1);
        
        if( rho(ii) < abs(eigen_2) )
            rho(ii) = abs(eigen_2);
        end
        
    end
    
    semilogx(dt_T, rho, 'LineWidth', 3);
    hold on;
end

grid on;

