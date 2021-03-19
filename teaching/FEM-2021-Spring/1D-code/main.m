clear all; clc;

% material properties and input data
kappa = 5.0;

% exact solution
exact = @(x) exp(x);

f = @(x) -1.0 * kappa * exp(x);
h = @(x) -kappa;
g = @(x) exp(1);

% domain
omega_l = 0.0;
omega_r = 1.0;

% interpolation degree
pp = 1;

% quadrature rule
nqp = 2;
[qp, wq] = Gauss(nqp, -1, 1);

% finite element mesh
nElem = 5; 

n_np = nElem + 1;
n_en = pp + 1;

IEN = zeros(n_en, nElem);

for ee = 1 : nElem
    for aa = 1 : n_en
        IEN(aa, ee) = ee - 1 + aa;
    end
end

% mesh is assumbed to have uniform size hh
hh = (omega_r - omega_l) / nElem; 

x_coor = omega_l : hh : omega_r;

% setup ID array based on the boundary condition
ID = 1 : n_np;
ID(end) = 0;

% Setup the stiffness matrix and load vector
% number of equations equals the number of nodes minus the number of
% Dirichlet nodes
n_eq = n_np - 1;

K = sparse(n_eq, n_eq);
F = zeros(n_eq, 1);

% Assembly the siffness matrix and load vector
for ee = 1 : nElem
    k_ele = zeros(n_en, n_en);
    f_ele = zeros(n_en, 1);
    
    x_ele = zeros(n_en, 1);
    for aa = 1 : n_en
        x_ele(aa) = x_coor( IEN(aa,ee) );
    end
    
    for qua = 1 : nqp
        x_qua = 0.0;
        dx_dxi = 0.0;
        for aa = 1 : n_en
            x_qua = x_qua + x_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
            dx_dxi = dx_dxi + x_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
        end
        dxi_dx = 1.0 / dx_dxi;
        
        for aa = 1 : n_en
            Na    = PolyBasis(pp, aa, 0, qp(qua));
            Na_xi = PolyBasis(pp, aa, 1, qp(qua));
            f_ele(aa) = f_ele(aa) + wq(qua) * Na * f(x_qua) * dx_dxi;
            for bb = 1 : n_en
                Nb_xi = PolyBasis(pp, bb, 1, qp(qua));
                k_ele(aa,bb) = k_ele(aa,bb) + wq(qua) * Na_xi * kappa * Nb_xi * dxi_dx;
            end
        end        
    end
    % end of the quadrature loop
    
    % distribute the entries to the global stiffness matrix and global load vector
    for aa = 1 : n_en
        LM_a = ID( IEN(aa, ee) );
        if LM_a > 0
            F(LM_a) = F(LM_a) + f_ele(aa);
            for bb = 1 : n_en
                LM_b = ID( IEN(bb, ee) );
                if LM_b > 0
                    K(LM_a, LM_b) = K(LM_a, LM_b) + k_ele(aa, bb);
                else
                    x_qua = x_coor( IEN(bb,ee) ); % obtain the Dirichlet node's physical coordinates
                    g_qua = g( x_qua ); % Obtain the boundary data at this point
                    F( LM_a ) = F( LM_a ) - k_ele(aa, bb) * g_qua;
                end
            end
        end
    end
    
    % Modify the load vector by the Natural BC
    % Note: for multi-dimensional cases, one needs to perform line or
    % surface integration for the natural BC.
    if ee == 1
        F( ID(IEN(1, ee)) ) = F( ID(IEN(1, ee)) ) + h( x_coor(IEN(1,ee)));
    end   
end

% Solve the stiffness matrix problem
disp = K \ F;

% Append the displacement vector by the Dirichlet data
disp = [ disp; g(omega_r) ];

% EOF