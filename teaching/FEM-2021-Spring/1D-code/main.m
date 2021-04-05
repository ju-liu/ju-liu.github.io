% clean the memory and the screen
clear all; clc;

% material properties and input data
kappa = 5.0;

% exact solution
exact = @(x) sin(x);
exact_x = @(x) cos(x);

f = @(x) kappa * sin(x);
h = @(x) -kappa;
g = @(x) sin(1);

% domain
omega_l = 0.0;
omega_r = 1.0;

% interpolation degree
pp = 5;

% number of elements
nElem = 4;

% quadrature rule
nqp = 4;
[qp, wq] = Gauss(nqp, -1, 1);

n_np = nElem * pp + 1; % number of nodal points
n_en = pp + 1;         % number of element nodes

IEN = zeros(n_en, nElem);

for ee = 1 : nElem
    for aa = 1 : n_en
        IEN(aa, ee) = (ee - 1) * pp + aa;
    end
end

% mesh is assumbed to have uniform size hh
hh = (omega_r - omega_l) / nElem;

x_coor = omega_l : (hh/pp) : omega_r;

% setup ID array based on the boundary condition
ID = 1 : n_np;
ID(end) = 0;

% Setup the stiffness matrix and load vector
% number of equations equals the number of nodes minus the number of
% Dirichlet nodes
n_eq = n_np - 1;

% Allocate an empty stiffness matrix and load vector
K = sparse(n_eq, n_eq);
F = zeros(n_eq, 1);

% Assembly the siffness matrix and load vector
for ee = 1 : nElem
    % Allocate zero element stiffness matrix and element load vector
    k_ele = zeros(n_en, n_en);
    f_ele = zeros(n_en, 1);
    
    x_ele = zeros(n_en, 1);
    for aa = 1 : n_en
        x_ele(aa) = x_coor( IEN(aa,ee) );
    end
    
    for qua = 1 : nqp
        % geometrical mapping
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
uh = K \ F;

% Append the displacement vector by the Dirichlet data
uh = [ uh; g(omega_r) ];

% Now we do the postprocessing
nqp = 6;
[qp, wq] = Gauss(nqp, -1, 1);

top = 0.0; bot = 0.0;
for ee = 1 : nElem
    for qua = 1 : nqp
        x_ele = zeros(n_en, 1);
        u_ele = zeros(n_en, 1);
        for aa = 1 : n_en
            x_ele(aa) = x_coor(IEN(aa, ee));
            u_ele(aa) = uh(IEN(aa, ee));
        end
        
        x = 0.0; dx_dxi = 0.0; duh_dxi = 0.0;
        for aa = 1 : n_en
            x = x + x_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
            dx_dxi = dx_dxi + x_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
            duh_dxi = duh_dxi + u_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
        end
        
        dxi_dx = 1.0 / dx_dxi;
        
        top = top + wq(qua) * (duh_dxi * dxi_dx - exact_x(x))^2 * dx_dxi;
        bot = bot + wq(qua) * exact_x(x)^2 * dx_dxi;
    end
end

top = sqrt(top);
bot = sqrt(bot);

error = top / bot;

% Remark:
%
% h-Fem or h-version fo the finite element method refers to maintaining the
% polynomial degree (probably low-order polynomial) and reducing the 
% element size. We note here h usually refers to the characteristic element
% size.
%
% p-FEM or p-version of the finite element method refers to a
% discretization strategy by elevating the polynomial degree with mesh
% fixed.
%
% hp-FEM referes to using a suitable combination of h-refinements and
% p-refinements, and it may give exponential convergence rates.
%
% EOF