function RdR = NURBSbasis(I,xi,degree,knots,w,n,computeMixedDerivs)
% This routine compute the nonzero NURBS functions 
% and corresponding derivatives at pt
% This routine can be generalized based on Algorithm A4.2 in The NURBS Book Piegl1997tnb available at
% https://books.google.ch/books?hl=en&lr=&id=7dqY5dyAwWkC&oi=fnd&pg=PA1&dq=the+nurbs+book+piegl+1997&ots=SACXTVU4qR&sig=skARkgLc0V95fQfCcQRm3qR9hwI#v=onepage&q=the%20nurbs%20book%20piegl%201997&f=false

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
if nargin < 7
    computeMixedDerivs = false;
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
for j = 1:d_p
    B{j} = BsplineBasis(I(:,j), xi(:,j), degree(j), knots{j}, n);
    if n > 0
        dW{j} = zeros(noxi,1,n,prec);
        dR{j} = zeros(noxi,n_en,n,prec);
    end
end

BB = ones([noxi,order],prec);
if n > 0
    dBB = cell(n,d_p);
end
for j = 1:d_p
    temp_i = ones(1,d_p);
    temp_i(j) = order(j);
    BB = BB.*reshape(B{j}(:,:,1),[noxi,temp_i]);
    if n > 0
        for i = 1:d_p
            for k = 1:n
                if j == 1   
                    dBB{k,i} = ones(noxi,1,prec);
                end
                if j == i
                    der_idx = k+1;
                else
                    der_idx = 1;
                end
                dBB{k,i} = dBB{k,i}.*reshape(B{j}(:,:,der_idx),[noxi,temp_i]);
            end
        end
    end
end

w = reshape(w,[numel(w)/n_en,order]); % Reshape and take into account multiple knot spans
W = sum(BB.*w,2:d_p+1);
if n > 0
    for j = 1:d_p
        for k = 1:n
            dW{j}(:,:,k) = sum(dBB{k,j}.*w,2:d_p+1);
        end
    end
end

R = BB.*w./W;
if n > 0
    dRtemp = cell(n+1,d_p); % non-mixed derivatives
    for k = 0:n
        for j = 1:d_p
            if k == 0
                dRtemp{k+1,j} = R;
            else
                dRtemp{k+1,j} = zeros([noxi,order],prec);
            end
        end
    end
    for k = 1:n
        for j = 1:d_p
            for i = 1:k
                dRtemp{k+1,j} = dRtemp{k+1,j} + nchoosek(k,i)*dW{j}(:,:,i).*dRtemp{k-i+1,j};
            end
            dRtemp{k+1,j} = (w.*dBB{k,j} - dRtemp{k+1,j})./W; % Generalization of Eq. (4.8) w.r.t. the number of parametric directions d_p (Eq. (4.8) has d_p=1)
        end
    end
end
if n > 1 && computeMixedDerivs
    noMixedDerivs = nchoosek(d_p,n);
    if n > 2
        warning('Mixed derivatives not implemented for order larger than 2')
    end
    d2BB = cell(1,noMixedDerivs);
    d2A = cell(1,noMixedDerivs);
    d2W = cell(1,noMixedDerivs);
    for i = 1:noMixedDerivs
        d2BB{i} = ones(noxi,1,prec);
        for j = 1:d_p
            temp_i = ones(1,d_p);
            temp_i(j) = order(j);
            if d_p == 3
                if j == i % Use ordering inspired by Voigt notation for stress tensors (11,22,33,23,13,12)
                    der_idx = 1;
                else
                    der_idx = 2;
                end
            elseif d_p == 2
                der_idx = 2;
            end
            d2BB{i} = d2BB{i}.*reshape(B{j}(:,:,der_idx),[noxi,temp_i]);
        end
        d2A{i} = w.*d2BB{i};
        d2W{i} = sum(d2BB{i}.*w,2:d_p+1);
    end

    % n = 2:
    % d_p = 2: {(1,2)}
    % d_p = 3: {(2,3), (1,3), (1,2)}
    % d_p = 4: {(1,4), (2,4), (3,4), (1,3), (2,3), (1,2)}
    dmRtemp = cell(1,noMixedDerivs);
    switch d_p
        case 2
            dmRtemp{1} = (d2A{1} - dW{1}(:,:,1).*dRtemp{2,2} - dW{2}(:,:,1).*dRtemp{2,1} - d2W{1}.*R)./W;
        case 3
            dmRtemp{1} = (d2A{1} - dW{2}(:,:,1).*dRtemp{2,3} - dW{3}(:,:,1).*dRtemp{2,2} - d2W{1}.*R)./W;
            dmRtemp{2} = (d2A{2} - dW{1}(:,:,1).*dRtemp{2,3} - dW{3}(:,:,1).*dRtemp{2,1} - d2W{2}.*R)./W;
            dmRtemp{3} = (d2A{3} - dW{1}(:,:,1).*dRtemp{2,2} - dW{2}(:,:,1).*dRtemp{2,1} - d2W{3}.*R)./W;
        otherwise
            error('Not implemented')
    end
    dmR = cell(1,noMixedDerivs); % Mixed derivatives only for n >= 2
    for j = 1:noMixedDerivs
        dmR{j} = reshape(dmRtemp{j},noxi,n_en);
    end
end

R = reshape(R,noxi,n_en);
if n > 0
    for j = 1:d_p
        for k = 1:n
            dR{j}(:,:,k) = reshape(dRtemp{k+1,j},noxi,n_en);
        end
    end
end

if n == 0
    RdR = {R};
elseif n == 1 || ~computeMixedDerivs
    RdR = [{R},dR];
else
    RdR = [{R},dR,dmR];
end



