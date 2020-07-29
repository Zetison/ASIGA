function I = radialIntegral(n, r_a, k, f, dimension, infiniteElementFormulation, type)


if dimension == 3
    switch infiniteElementFormulation
        case {'PGC', 'BGC',}
            if type == 1
                % I = integral(1/r^n) from r_a to infinity
                I = 1/r_a^(n-1)/(n-1);
            elseif type == 2
                % I = integral(1/(r^n*(r^2-f^2)) from r_a to infinity
                if f == 0
                    I = 1/r_a^(n+1)/(n+1);
                else
                    I = 0;
                    I_j = inf;
                    j = 0;
                    fDr_a2 = (f/r_a)^2;
                    pFact = 1;
                    while abs(I_j/I) > eps
                        I_j = pFact/(2*j+n+1);
                        I = I + I_j;
                        j = j + 1;
                        pFact = pFact*fDr_a2;
                    end
                    I = I/r_a^(n+1);
                end
            end
        case {'PGU', 'BGU'}
            if type == 1
                % I = integral(exp(2*1i*k*r)/r^n) from r_a to infinity
                I = 1/r_a^(n-1)*expint_n(-2*1i*k*r_a,n);
            elseif type == 2
                % I = integral(exp(2*1i*k*r)/(r^n*(r^2-f^2)) from r_a to infinity
                if f == 0
                    I = 1/r_a^(n+1)*expint_n(-2*1i*k*r_a,n+2);
                else
                    I = 0;
                    I_j = inf;
                    j = 0;
                    fDr_a2 = (f/r_a)^2;
                    pFact = 1;
                    while abs(I_j/I) > eps
                        I_j = pFact*expint_n(-2*1i*k*r_a,2*j+n+2);
                        I = I + I_j;
                        j = j + 1;
                        pFact = pFact*fDr_a2;
                    end
                    I = I/r_a^(n+1);
                end
            end
    end
elseif dimension == 2
    switch infiniteElementFormulation
        case {'PGC', 'BGC',}
            if type == 1
                % I = integral(1/(r^n*sqrt(r^2-f^2)) from r_a to infinity
                if f == 0
                    I = 1/r_a^n/n;
                else
                    I = 0;
                    I_j = inf;
                    j = 0;
                    fDr_a2 = (f/r_a)^2;
                    pFact = 1;
                    a_j = 1;
                    while abs(I_j/I) > eps
                        I_j = pFact*a_j/(2*j+n);
                        I = I + I_j;
                        j = j + 1;
                        pFact = pFact*fDr_a2;
                        a_j = a_j/2*(2*j-1)/j;
                    end
                    I = I/r_a^n;
                end
            elseif type == 2
                I = NaN; % This type does not exist
            end
        case {'PGU', 'BGU'}
            if type == 1
                % I = integral(exp(2*1i*k*r)/(r^n*sqrt(r^2-f^2)) from r_a to infinity
                if f == 0
                    I = 1/r_a^n*expint_n(-2*1i*k*r_a,n+1);
                else
                    I = 0;
                    I_j = inf;
                    j = 0;
                    fDr_a2 = (f/r_a)^2;
                    pFact = 1;
                    a_j = 1;
                    while abs(I_j/I) > eps
                        I_j = pFact*a_j*expint_n(-2*1i*k*r_a,2*j+n+1);
                        I = I + I_j;
                        j = j + 1;
                        pFact = pFact*fDr_a2;
                        a_j = a_j/2*(2*j-1)/j;
                    end
                    I = I/r_a^n;
                end
            elseif type == 2
                % I = integral(exp(2*1i*k*r)/((r^2-f^2)^(3/2)) from r_a to infinity
                if f == 0
                    I = 1/r_a^2*expint_n(-2*1i*k*r_a,3);
                else
                    I = 0;
                    I_j = inf;
                    j = 0;
                    fDr_a2 = (f/r_a)^2;
                    pFact = 1;
                    a_j = 1;
                    while abs(I_j/I) > eps
                        I_j = pFact*a_j*(2*j+1)*expint_n(-2*1i*k*r_a,2*j+3);
                        I = I + I_j;
                        j = j + 1;
                        pFact = pFact*fDr_a2;
                        a_j = a_j/2*(2*j-1)/j;
                    end
                    I = I/r_a^2;
                end            
            end
    end
end

