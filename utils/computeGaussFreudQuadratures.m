function [Q,W] = computeGaussFreudQuadratures(n,r)
useSymbolic = isa(r,'sym');
if useSymbolic
    Eps = 10^(-digits);
else
    Eps = 10*eps;
end
maxItrs = 1000;

Q = cell(n+1,1);
W = cell(n,1);
A = zeros(n,1);
gamma = zeros(n,1);
if useSymbolic
    gamma = vpa(gamma);
    A = vpa(A);
end

[Bprev, g] = getFreudCoeffs(0, r);
for m = 1:n
    m
    [B, gamma(m)] = getFreudCoeffs(m, r);
    if r >= 4
        regF = @(n) 0.948+(-0.1461+0.008928*n)*(exp(-0.01558*(n-4)^2)*heaviside(n-4)+1-heaviside(n-4));
    else
        regF = @(n) 0.9347+(-0.1461+0.008928*n)*(exp(-0.01558*(n-4)^2)*heaviside(n-4)+1-heaviside(n-4));
    end
    if m == 2
        bndrs = [0, Q{m-1}, 3*Q{m-1}(m-1)];
    elseif m > 2
        bndrs = [0, Q{m-1}, 2*Q{m-1}(m-1)-Q{m-1}(m-2)];
    end
    for i = 1:m
        for j = 1:4
            if m == 1
                x_0 = 0;
            else
                if j == 1
                    x_0 = bndrs(i+1)*regF(m-1);
                elseif j == 2
                    x_0 = 1/4*bndrs(i)+3/4*bndrs(i+1);   
                elseif j == 3
                    x_0 = 3/4*bndrs(i)+1/4*bndrs(i+1);    
                elseif j == 4
                    x_0 = 1/2*bndrs(i)+1/2*bndrs(i+1);                
                end
            end
            Q{m}(i) = newtonsMethod(@(x)P_(B,x,m),@(x)dPdx_(B,x,m),x_0,maxItrs,10^(-log10(Eps)/4)*Eps);
            if ~(m == 1) && bndrs(i) < Q{m}(i) && Q{m}(i) < bndrs(i+1) 
                break
            end
        end
    end
    [~, I] = uniquetol(double(Q{m}),10*eps);
    Q{m} = Q{m}(I);
    if length(Q{m}) < m
        warning(['Not all quadrature points were found when m = ' num2str(m) '.'])
    end
    A(m) = (-1)^m/prod(Q{m});
    if m == 1
        a_m = A(m);
    else
        a_m = A(m)/A(m-1);
        g = gamma(m-1);
    end
    W{m} = zeros(1,m);
    if useSymbolic
        W{m} = vpa(W{m});
    end
    for i = 1:m  
        x_i = Q{m}(i);          
        W{m}(i) = a_m*g/(dPdx_(B,x_i,m)*P_(Bprev,x_i,m-1));
    end
    Bprev = B;
end

function P = P_(B,x,n)
P = 0;
for m = 0:n
    P = P + B(m+1)*x.^m;
end

function dPdx = dPdx_(B,x,n)
dPdx = 0;
for m = 1:n
    dPdx = dPdx + m*B(m+1)*x.^(m-1);
end



