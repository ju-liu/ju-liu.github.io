% =========================================================================
% This is the shape function routine for one-dimensional finite element code
%
% In this function, we give the Lagrange type basis function over the
% computational domain [-1,1]. Interior nodes are uniformly distributed in
% the reference domain.
%
% degree: the degree of our interpolation basis function space. In our code
%         we give the basis function up to degree six.
% i     : the number of the basis function. i takes value from 1 to degree+1.
% der   : if der == 0, return the value of the basis function;
%         if der == 1, return the 1st derivative of the basis function;
%         if der == 2, return the 2nd derivative of the basis funciton.
% x     : the point we which to perform the evaluation.
%
% Output: the value of the basis function or the 1st and 2nd derivatives of 
%         the basis function.
% -------------------------------------------------------------------------
% By Ju Liu, CAM student, 2009 Dec. 24th.
% =========================================================================

function poly = PolyBasis(degree , i , der , x)
switch degree
    % linear basis function
    case 1
        if i == 1
            if der == 0
                poly = 0.5 * (1-x);
            elseif der == 1
                poly = -0.5;
            elseif der == 2
                poly = 0;
            end
        end
        if i == 2
            if der == 0
                poly = 0.5 * (1+x);
            elseif der == 1
                poly = 0.5;
            elseif der == 2
                poly = 0;
            end
        end
        
        % quadratic basis function
    case 2
        if i == 1
            if der == 0
                poly = 0.5 * x * (x-1);
            elseif der == 1
                poly = x-0.5;
            elseif der == 2
                poly = 1;
            end
        end
        if i == 2
            if der == 0
                poly = 1 -x^2;
            elseif der == 1
                poly = -2 * x;
            elseif der == 2
                poly = -2;
            end
        end
        if i == 3
            if der == 0
                poly = 0.5 * x * (x+1);
            elseif der == 1
                poly = x+0.5;
            elseif der == 2
                poly = 1;
            end
        end
        
        % Cubic basis function
    case 3
        if i == 1
            if der == 0
                poly = -9*(x-(1/3)) * (x+(1/3)) * (x-1)/16;
            elseif der == 1
                poly = -9*(2*x*(x-1)+x^2-(1/9))/16;
            elseif der == 2
                poly = -27/8*x+9/8;
            end
        end
        if i == 2
            if der == 0
                poly = 27 *(x^2-1)*(x-(1/3))/16;
            elseif der == 1
                poly = 27 * (2*x*(x-(1/3))+x^2-1)/16;
            elseif der == 2
                poly = 81/8*x-9/8;
            end
        end
        if i == 3
            if der == 0
                poly = -27 * (x^2-1)*(x+(1/3))/16;
            elseif der == 1
                poly = -27 * (2*x*(x+(1/3))+x^2-1)/16;
            elseif der == 2
                poly = -81/8*x-9/8;
            end
        end
        if i == 4
            if der == 0
                poly = 9*(x+1)*(x^2-(1/9))/16;
            elseif der == 1
                poly = 9*(x^2-(1/9)+2*x*(x+1))/16;
            elseif der == 2
                poly = 27/8*x+9/8;
            end
        end
        
        % quartic basis function
    case 4
        if i == 1
            if der == 0
                poly = 2*x*(x^2-(1/4))*(x-1)/3;
            elseif der == 1
                poly = 2*((x^2-(1/4))*(x-1)+2*x^2*(x-1)+x*(x^2-(1/4)))/3;
            elseif der == 2
                poly = 4*x*(x-1)+4*x^2-1/3;
            end
        end
        if i == 2
            if der == 0
                poly = -8*x*(x^2-1)*(x-0.5)/3;
            elseif der == 1
                poly = -8*((x^2-1)*(x-0.5)+x^2*(2*x-1)+x*(x^2-1))/3;
            elseif der == 2
                poly = -16/3*x*(x-1/2)-16*x^2+16/3-16/3*x*(2*x-1);
            end
        end
        if i == 3
            if der == 0
                poly = 4*(x^2-1)*(x^2-0.25);
            elseif der == 1
                poly = 4*(2*x*(x^2-0.25)+2*x*(x^2-1));
            elseif der == 2
                poly = 48*x^2-10;
            end
        end
        if i == 4
            if der == 0
                poly = -8*x*(x^2-1)*(x+0.5)/3;
            elseif der == 1
                poly = -8*((x^2-1)*(x+0.5)+x^2*(2*x+1)+x*(x^2-1))/3;
            elseif der == 2
                poly = -16/3*x*(x+1/2)-16*x^2+16/3-16/3*x*(2*x+1);
            end
        end
        if i == 5
            if der == 0
                poly = 2*x*(x^2-0.25)*(x+1)/3;
            elseif der == 1
                poly = 2*((x^2-0.25)*(x+1)+2*x^2*(x+1)+x*(x^2-0.25))/3;
            elseif der == 2
                poly = 4*x*(x+1)+4*x^2-1/3;
            end
        end
        
        % quintic basis function
    case 5
        if i == 1
            if der == 0
                poly =-625*(x^2-(9/25))*(x^2-(1/25))*(x-1)/768;
            elseif der == 1
                poly = -3125/768*x^4+625/192*x^3+125/128*x^2-125/192*x-3/256;
            elseif der == 2
                poly = -3125/192*x^3+625/64*x^2+125/64*x-125/192;
            end
        end
        
        if i == 2
            if der == 0
                poly = 3125/768*(x+1)*(x+1/5)*(x-1/5)*(x-3/5)*(x-1);
            elseif der == 1
                poly = 15625/768*x^4-625/64*x^3-1625/128*x^2+325/64*x+125/768;
            elseif der == 2
                poly = 15625/192*x^3-1875/64*x^2-1625/64*x+325/64;
            end
        end
        
        if i == 3
            if der == 0
                poly = -3125/384*(x+1)*(x+3/5)*(x-1/5)*(x-3/5)*(x-1);
            elseif der == 1
                poly = -15625/384*x^4+625/96*x^3+2125/64*x^2-425/96*x-375/128;
            elseif der == 2
                poly = -15625/96*x^3+625/32*x^2+2125/32*x-425/96;
            end
        end
        
        if i == 4
            if der == 0
                poly = 3125/384*(x+1)*(x+3/5)*(x+1/5)*(x-3/5)*(x-1);
            elseif der == 1
                poly = 15625/384*x^4+625/96*x^3-2125/64*x^2-425/96*x+375/128;
            elseif der == 2
                poly = 15625/96*x^3+625/32*x^2-2125/32*x-425/96;
            end
        end
        
        if i == 5
            if der == 0
                poly = -3125/768*(x+1)*(x+3/5)*(x+1/5)*(x-1/5)*(x-1);
            elseif der == 1
                poly = -15625/768*x^4-625/64*x^3+1625/128*x^2+325/64*x-125/768;
            elseif der == 2
                poly = -15625/192*x^3-1875/64*x^2+1625/64*x+325/64;
            end
        end
        
        if i == 6
            if der == 0
                poly = 625/768*(x+1)*(x+3/5)*(x+1/5)*(x-1/5)*(x-3/5);
            elseif der == 1
                poly = 3125/768*x^4-125/128*x^2+3/256+625/192*x^3-125/192*x;
            elseif der == 2
                poly = 3125/192*x^3-125/64*x+625/64*x^2-125/192;
            end
        end
        
    case 6
        if i == 1
            if der == 0
                poly = 81/80*(x+2/3)*(x+1/3)*x*(x-1/3)*(x-2/3)*(x-1);
            elseif der == 1
                poly = 1/80*(6*x-1)*(81*x^4-54*x^3-39*x^2+16*x+4);
            elseif der == 2
                poly = 243/40*x^4-81/20*x^3-117/40*x^2+6/5*x+3/10+(3/40*x-1/80)*(324*x^3-162*x^2-78*x+16);
            end
        end
        
        if i == 2
            if der == 0
                poly = -243/40*(x+1)*(x+1/3)*x*(x-1/3)*(x-2/3)*(x-1);
            elseif der == 1
                poly = -27/20*x-27/2*x^2+81/4*x^4-729/20*x^5+27*x^3+9/20;
            elseif der == 2
                poly = -27/20-27*x+81*x^3-729/4*x^4+81*x^2;
            end
        end
        
        if i == 3
            if der == 0
                poly = 243/16*(x+1)*(x+2/3)*x*(x-1/3)*(x-2/3)*(x-1);
            elseif der == 1
                poly = 27/2*x+351/16*x^2-405/16*x^4-9/4+729/8*x^5-351/4*x^3;
            elseif der == 2
                poly = 27/2+351/8*x-405/4*x^3+3645/8*x^4-1053/4*x^2;
            end
        end
        
        if i == 4
            if der == 0
                poly = -81/4*(x+1)*(x+1/3)*(x+2/3)*(x-1/3)*(x-2/3)*(x-1);
            elseif der == 1
                poly = -49/2*x-243/2*x^5+126*x^3;
            elseif der == 2
                poly = -49/2-1215/2*x^4+378*x^2;
            end
        end
        
        if i == 5
            if der == 0
                poly = 243/16*(x+1)*(x+1/3)*(x+2/3)*x*(x-2/3)*(x-1);
            elseif der == 1
                poly = 27/2*x-351/16*x^2+405/16*x^4+729/8*x^5+9/4-351/4*x^3;
            elseif der == 2
                poly = 27/2-351/8*x+405/4*x^3+3645/8*x^4-1053/4*x^2;
            end
        end
        
        if i == 6
            if der == 0
                poly = -243/40*(x+1)*(x+1/3)*(x+2/3)*x*(x-1/3)*(x-1);
            elseif der == 1
                poly = -27/20*x+27/2*x^2-81/4*x^4-729/20*x^5+27*x^3-9/20;
            elseif der == 2
                poly = -27/20+27*x-81*x^3-729/4*x^4+81*x^2;
            end
        end
        
        if i == 7
            if der == 0
                poly = 81/80*(x+1)*(x+1/3)*(x+2/3)*x*(x-1/3)*(x-2/3);
            elseif der == 1
                poly = 1/80*(6*x+1)*(81*x^4+54*x^3-39*x^2-16*x+4);
            elseif der == 2
                poly = 243/40*x^4+81/20*x^3-117/40*x^2-6/5*x+3/10+(3/40*x+1/80)*(324*x^3+162*x^2-78*x-16);
            end
        end
end