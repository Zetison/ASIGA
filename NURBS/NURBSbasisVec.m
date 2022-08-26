function RdR = NURBSbasisVec(I,xi,degree,knots,w,n)
% This routine compute the nonzero NURBS functions 
% and corresponding derivatives at pt
% This routine can be generalized based on Algorithm A4.2 in The NURBS Book Piegl1997tnb available at
% https://doi.org/10.1016/S0167-8396(98)00028-4

% Input
%       xi:                 evaluation points
%       degree:             NURBS degrees
%       knots:              knot vectors
%       w:                  NURBS weights

% Output
%       R, dR{1}, dR{2}, dR{3}:      array of the (p+1)(q+1)(r+1) NURBS 
%                                       functions and its derivatives which 
%                                       are nonzero at (xi,eta,zeta)
if nargin < 6
    n = 1;
end
d_p = numel(degree);
prec = class(xi);
noxi = size(xi,1);
order = degree+1;
n_en = prod(order);
dd = ones(d_p);
B = cell(1,d_p);
if n > 0
    dW = cell(1,d_p);
    dR = cell(1,d_p);
else
    dR = {};
end
for i = 1:d_p
    dd(i,i) = 2;
    B{i} = BsplineBasisVec(I(:,i), xi(:,i), degree(i), knots{i}, n);
    if n > 0
        dW{i} = zeros(noxi,1,n,prec);
        dR{i} = zeros(noxi,n_en,n,prec);
    end
end

if n > 1
    error('Not implemented yet')
end
BB = ones(size(w),prec);
if n > 0
    dBB = cell(1,d_p);
end
for i = 1:d_p
    temp_i = ones(1,d_p);
    temp_i(i) = degree(i)+1;
    BB = BB.*reshape(B{i}(:,:,1),[noxi,temp_i]);
    if n > 0
        for j = 1:d_p
            if i == 1   
                dBB{j} = ones(size(w),prec);
            end
            dBB{j} = dBB{j}.*reshape(B{i}(:,:,dd(j,i)),[noxi,temp_i]);
        end
    end
end
W = sum(BB.*w,2:d_p+1);
if n > 0
    for i = 1:d_p
        for j = 1:n
            dW{i}(:,:,j) = sum(dBB{i}.*w,2:d_p+1);
        end
    end
end
fact = w./W;
R = reshape(BB.*fact,noxi,n_en);
if n > 0
    for j = 1:d_p
        dR{j}(:,:,1) = reshape((dBB{j} - BB.*dW{j}(:,:,1)./W).*fact,noxi,n_en);
    end
end
RdR = [{R},dR];