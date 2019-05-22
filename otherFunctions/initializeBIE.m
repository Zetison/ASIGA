function [constants, integrals] = initializeBIE(psiType,useRegul,x,nx,k,model)
switch psiType
    case 1
        switch model
            case {'Cube','Cube_P'}
                x1 = zeros(size(x));
            otherwise
                x1 = x-nx;
        end
        C1 = norm(x-x1);
        C2 = dot(x-x1, nx);
        constants = {x1 C1 C2};
    case 2
        switch model
            case {'Cube','Cube_P'}
                x1 = zeros(size(x));
                x2 = 0.5*x;
            otherwise
                x1 = x - 0.5*nx;
                x2 = x - nx;
        end
        r1x = norm(x1-x);
        r2x = norm(x2-x);
        C2 = (1i*k*r2x-1)/r2x^2*dot(x2-x,nx) - (1i*k*r1x-1)/r1x^2*dot(x1-x,nx);
        C1 = 1 - r2x^2*(1i*k*r1x-1)*dot(x1-x,nx)/(r1x^2*(1i*k*r2x-1)*dot(x2-x,nx));

        if abs(C2) < 1e-4 || abs(C1) < 1e-4
            error('Choose x1 and x2 more visely')
        end
        constants = {x1 x2 C1 C2 Phi_k(r1x,k) Phi_k(r2x,k)};
    case 3
        if abs(nx(1)) < 1/sqrt(2)
            d1 = sqrt(3)/2*cross([1,0,0],nx)/sqrt(1-nx(1)^2) - nx/2;
        else
            d1 = sqrt(3)/2*cross([0,1,0],nx)/sqrt(1-nx(2)^2) - nx/2;
        end
        d2 = d1+nx;
        constants = {d1 d2};
    otherwise
        constants = NaN;
end
integrals = cell(4,1);
if useRegul
    integrals{1} = complex(0);	% Psi1_integral
    integrals{2} = complex(0);	% Psi2_integral
    integrals{3} = complex(0);	% dPsi1dny_integral
    integrals{4} = complex(0);	% dPsi2dny_integral
else
    integrals{1} = complex(0);
    integrals{2} = complex(0);
    integrals{3} = complex(zeros(3,1));
end