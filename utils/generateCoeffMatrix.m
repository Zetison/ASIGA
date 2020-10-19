function [D,Dt,y_n] = generateCoeffMatrix(varCol)

r_a = varCol.r_a;

if varCol.IElocSup
    p_ie = varCol.p_ie;
    N = p_ie;
    if strcmp(varCol.IEbasis,'Lagrange')
        n_n = varCol.N;
        kappa = 1:n_n;
        if varCol.s_ie == 1
            z = 1 + (1-kappa)/n_n;
        elseif varCol.s_ie == 2
            z = 1 + kappa.*(1-kappa)/(n_n*(n_n+1));
        else
            warning('Not implemented, using s = 2 instead')
            z = 1 + (1-kappa)/n_n;
        end
        y_n = 1./z;
        x = z(end-p_ie+1:end)/z(end-p_ie+1);
    else
        y_n = NaN;
    end
else
    N = varCol.N;
%     a = 0.65;
%     b = 9;
%     logb = @(x,b) log(x)/log(b);
    m = 1:N;
    x = 1./m;
%     x(m) = 1./(1+(1-m)/N);
%     x(m) = (N-m+1)/N;
%     x(m) = a+(1-a)*(N-m+1)/N;
%     x(m) = 1./log10(10*m);
%     x(m) = 1./10.^(m-1);
%     x(m) = 1./logb(b*m,b);
%     GLL = gaussLobattoLegendreQuad(N+1);
%     x = 1./(0.5-0.5*GLL(1:end-1));
    y_n = NaN;
end
switch varCol.IEbasis
    case 'Standard'
        error('It is assumed that the radial basis functions satisfies the Kronecker delta property at the artificial boundary')
        k = varCol.k;
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
        B = zeros(N);
        Bt = zeros(N);
        E = zeros(N);
        k = varCol.k;
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
            if varCol.IElocSup
                E(m,m) = 1;
            else
                E(m,m) = exp(1i*k*r_a*(1-1/x(m)));
            end
        end
        D = E/B;
        Dt = E/Bt;   
%         Binv = invVandermonde(x.');   
%         D = E*Binv;
%         Dt = E*Binv;   
end