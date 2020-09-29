function [D,Dt] = generateCoeffMatrix(varCol)

N = varCol.N;
switch varCol.IEbasis
    case 'Standard'
        error('It is assumed that the radial basis functions satisfies the Kronecker delta property at the artificial boundary')
        k = varCol.k;
        r_a = varCol.r_a;
        D = eye(N)*exp(1i*k*r_a);
        Dt = D;
    case 'Chebyshev'        
        D = zeros(N);
        for i = 1:N
            D(i,1) = (-1)^(i+1);
        end
        if N > 1
            D(2,2) = 2;
            for i = 3:N
                for m = 2:i-1
                    D(i,m) = 4*D(i-1,m-1) - 2*D(i-1,m) - D(i-2,m);
                end
                D(i,i) = 4*D(i-1,i-1);
            end
        end
        
        D(2:end,1) = D(2:end,1) - 1;
        Dt = D;
    case 'Bernstein'        
        D = zeros(N);
        p = N-1;
        for i = 0:p
            for j = 0:p-i
                D(i+1,i+1+j) = (-1)^j*nchoosek(p,i)*nchoosek(p-i,j);
            end
        end
        D = flipud(D);
        Dt = D;
    case 'Lagrange'
        x = zeros(N,1);
%         a = 0.65;
%         b = 9;
%         logb = @(x,b) log(x)/log(b);
        for m = 1:N
            x(m) = 1/m;
%             x(m) = (N-m+1)/N;
%             x(m) = a+(1-a)*(N-m+1)/N;
%             x(m) = 1/log10(10*m);
%             x(m) = 1/10^(m-1);
%             x(m) = 1/logb(b*m,b);
        end
        B = zeros(N);
        Bt = zeros(N);
        E = zeros(N);
        k = varCol.k;
        r_a = varCol.r_a;
        for m = 1:N
            for l = 1:N
                B(m,l) = x(l)^m;
                switch varCol.formulation
                    case {'PGC', 'PGU'}
                        Bt(m,l) = x(l)^(m+2);
                    case {'BGC', 'BGU'}
                        Bt(m,l) = B(m,l);
                end
            end
            E(m,m) = exp(1i*k*r_a*(1-1/x(m)));
        end
        D = E/B;
        Dt = E/Bt;      
end