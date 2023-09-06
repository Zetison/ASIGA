function RdR = NURBSbasis(I,xi,degree,knots,w,n)
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
B = cell(1,d_p);
if n > 0
    dW = cell(1,d_p);
    dR = cell(1,d_p);
end
for i = 1:d_p
    B{i} = BsplineBasisVec(I(:,i), xi(:,i), degree(i), knots{i}, n);
    if n > 0
        dW{i} = zeros(noxi,1,n,prec);
        dR{i} = zeros(noxi,n_en,n,prec);
    end
end

BB = ones([noxi,order],prec);
if n > 0
    dBB = cell(n,d_p);
end
for i = 1:d_p
    temp_i = ones(1,d_p);
    temp_i(i) = order(i);
    BB = BB.*reshape(B{i}(:,:,1),[noxi,temp_i]);
    if n > 0
        for j = 1:d_p
            for k = 1:n
                if i == 1   
                    dBB{k,j} = ones(noxi,1,prec);
                end
                if i == j
                    der_idx = k+1;
                else
                    der_idx = 1;
                end
                dBB{k,j} = dBB{k,j}.*reshape(B{i}(:,:,der_idx),[noxi,temp_i]);
            end
        end
    end
end

w = reshape(w,[numel(w)/n_en,order]); % Reshape and take into account multiple knot spans
W = sum(BB.*w,2:d_p+1);
if n > 0
    for i = 1:d_p
        for k = 1:n
            dW{i}(:,:,k) = sum(dBB{k,i}.*w,2:d_p+1);
        end
    end
end

R = BB.*w./W;
if n > 0
    dRtemp = cell(n,d_p);
    for k = 0:n
        for i = 1:d_p
            if k == 0
                dRtemp{k+1,i} = R;
            else
                dRtemp{k+1,i} = zeros([noxi,order],prec);
            end
        end
    end
    for k = 1:n
        for i = 1:d_p
            for j = 1:k
                dRtemp{k+1,i} = dRtemp{k+1,i} + nchoosek(k,j)*dW{i}(:,:,j).*dRtemp{k-j+1,i};
            end
            dRtemp{k+1,i} = (w.*dBB{k,i} - dRtemp{k+1,i})./W;
        end
    end
    for i = 1:d_p
        for k = 1:n
            dR{i}(:,:,k) = reshape(dRtemp{k+1,i},noxi,n_en);
        end
    end
end

R = reshape(R,noxi,n_en);

if n == 0
    RdR = {R};
else
    RdR = [{R},dR];
end



