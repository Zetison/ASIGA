function B = bernsteinBasis(x,p,D,scaling)
if nargin < 4
    scaling = 1;
end
if D == 0
    B = zeros(numel(x),p+1,class(x));
else
    B = zeros(numel(x),p+1,D,class(x));
end

pscl_temp = pascal(p+1,class(x));
pscl = cell(p+1,1);
for i = 1:p+1
    pscl{i} = diag(fliplr(pscl_temp),p+1-i);
end
fact = ones(D+1,1,class(x));
if D > 0
    fact(2) = p;
    for i = 3:D+1
        fact(i) = fact(i-1)*(p-i+2);
    end
end
for m = 0:D
    for i = 0:p
        temp = zeros(size(x),class(x));
        for k = max(0,i+m-p):min(i,m)
%             temp = temp + (-1)^(k+m)*nchoosek(m,k)*B_(i-k,p-m,x);
            temp = temp + (-1)^(k+m)*pscl{m+1}(k+1)*B_(i-k,p-m,x,pscl{p+1-m});
%             temp = temp + (-1)^(k+m)*factorial(m)/factorial(k)/factorial(m-k)*B_(i-k,p-m,x);
        end
%         B(:,i+1,m+1) = factorial(p)/factorial(p-m)*temp*scaling^m;
        B(:,i+1,m+1) = fact(m+1)*temp*scaling^m;
    end
end

function B = B_(i,p,x,pscl)

B = pscl(i+1)*(1-x).^(p-i).*x.^i;
% B = nchoosek(p,i)*(1-x).^(p-i).*x.^i;