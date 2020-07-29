function I = radialIntegral2(n, ra, k, Upsilon, infiniteElementFormulation, type)

switch infiniteElementFormulation
    case {'PGC', 'BGC',}
        if type == 1
%             I = integral(1/rho^n) from 1 to infinity
            I = 1/(n-1);
        elseif type == 2
%             I = integral(1/(rho^(n-1)*(rho^2-varrho^2)) from 1 to infinity
            if Upsilon == 0
                I = 1/n;
            else
                varrho = Upsilon./ra;
%                 I = hypergeom([1, n/2], (n+2)/2, varrho^2)/n;
                I = 0;
                I_j = inf;
                j = 0;
                pFact = ones(size(ra));
                varrho2 = varrho.^2;
                while any(abs(I_j./I) > eps)
                    I_j = pFact/(2*j+n);
                    I = I + I_j;
                    j = j + 1;
                    pFact = pFact.*varrho2;
                end
            end
        end
    case {'PGU', 'BGU'}
        if type == 1
%             I = integral(exp(-z*r)/r^n) from 1 to infinity
            z = -2*1i*k*ra;
            I = expint_n(z,n);
        elseif type == 2
%               I = integral(exp(-z*rho)/(rho^(n-1)*(rho^2-varrho^2)) from 1 to infinity
            if Upsilon == 0
                z = -2*1i*k*ra;
                I = expint_n(z,n+1);
            else
                z = -2*1i*k*ra;
                indices_1 = abs(z) > -inf; % indices for first series representation    
%                 indices_1 = abs(z) > 1; % indices for first series representation    
%                 indices_2 = abs(z) <= 1; % indices for second series representation
                I = zeros(size(z));
                if any(indices_1)
                    z = -2*1i*k*ra(indices_1);
                    varrho = Upsilon./ra(indices_1);
                    I(indices_1) = 0;
                    I_j = inf;
                    j = 0;
                    varrho2 = varrho.^2;
                    pFact = 1;
                    while any(abs(I_j./I(indices_1)) > eps)
                        I_j = pFact.*expint_n(z(indices_1),2*j+n+1);
                        I(indices_1) = I(indices_1) + I_j;
                        j = j + 1;
                        pFact = pFact.*varrho2;
                    end
                end
%                 if any(indices_2)
%                     z = -2*1i*k*ra(indices_2);
%                     varrho = Upsilon./ra(indices_2);
%                     gamma = 0.5772156649015328606065120900824;
%                     I(indices_2) = -1./(2*varrho.^n).*(exp(-varrho.*z).*(log((1-varrho).*z)+gamma) + (-1)^n*exp(varrho.*z).*(log((1+varrho).*z)+gamma));
%                     m = 0;
%                     pFact = 1;
%                     term = inf;
%                     while any(abs(term./I(indices_2)) > eps)
%                         D_m = 0;
%                         for j = 0:floor(n/2)-1
%                             if m == n-2*j-2
%                                 temp2 = log(z) + gamma - sum(1./(1:m));
%                             else
%                                 temp2 = 1/(m-n+2*j+2);
%                             end
%                             D_m = D_m + 1./varrho.^(2*(j+1)).*temp2;
%                         end
%                         if m ~= 0
%                             C_m = -(exp(-varrho.*z).*(1-varrho).^m+(-1)^n.*exp(varrho.*z).*(1+varrho).^m)./(2*m*varrho.^n);
%                         else
%                             C_m = 0;
%                         end
%                         term = pFact.*(C_m+D_m);
%                         I(indices_2) = I(indices_2) + term;
%                         m = m + 1;
%                         pFact = -z.*pFact/m;
%                     end
%                 end
            end
        end
end
