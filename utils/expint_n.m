function E_n = expint_n(x,n)
% Computes integral(exp(-x*y)/y^n dy, y=1..inf)

indices_s = abs(x) < 1; % indices for series representation    
indices_f = abs(x) >= 1; % indices for continued fraction representation    
E_n = zeros(size(x));


if any(indices_s)
    if n ~= 1
        E_n(indices_s) = 1/(n-1);
    else 
        E_n(indices_s) = -log(x(indices_s)) + psi(1);
    end
    m = 1;
    pFact = 1;
    term = inf;
    while any(abs(term./E_n(indices_s)) > eps)
        pFact = -x(indices_s)*pFact/m;
        if m ~= n-1
            term = -pFact/(m-n+1);
        else
            term = pFact*(-log(x(indices_s)) + psi(n));
        end
        E_n(indices_s) = E_n(indices_s) + term;
        m = m + 1;
    end
end
if any(indices_f)
    %     tiny = eps*mean(min(1,abs(x)))*1e-20;
    tiny = eps*mean(abs(x(indices_f)+n))*1e-20;
    f = tiny;
    C = f;
    D = 0;
    j = 1;
    Delta = zeros(size(x(indices_f)));
    while abs(Delta - 1) > eps
    %         if mod(j,2) == 1
    %             b = x;
    %             if j == 1
    %                 a = 1;
    %             else
    %                 a = (j-1)/2;
    %             end
    %         else
    %             b = 1;
    %             a = n + j/2 - 1;
    %         end
        if j == 1
            a = 1;
        else
            a = -(j-1)*(n+j-2);
        end
        b = x(indices_f) + n + 2*(j - 1);

        D = b + a*D;
        if D == 0
            D = tiny;
        end
        C = b+a./C;
        if C == 0
            C = tiny;
        end
        D = 1./D;
        Delta = C.*D;
        f = f.*Delta;
        j = j+1;
    end
    E_n(indices_f) = exp(-x(indices_f)).*f; 
end

function psi_n = psi(n)

gamma = 0.5772156649015328606065120900824;

psi_n = -gamma;
for m = 1:n-1
    psi_n = psi_n + 1/m;
end


