function Y = getInterpolatingHermite(x,X,n_D,plotPolynomials)
if nargin < 4
    plotPolynomials = false;
end
P = numel(x);
Y = zeros(P*n_D,numel(X),class(x));
ii = reshape(repmat(1:P,n_D,1),n_D*P,1);
nn = reshape(repmat(1:n_D,P,1).',n_D*P,1);
parfor counter = 1:P*n_D
% for counter = 1:P*n_D
    mp.Digits(400);
    i = ii(counter);
    n = nn(counter);
    A = [x, zeros(P,n_D,class(x))];
    A(i,n+1) = 1;
    [a, xj] = difftable(A);
    if 0
        YY = zeros(numel(a),numel(X),class(x));
        YY(1,:) = a(1)*ones(size(X,1),size(X,2),class(x));
        for m = 2:numel(a)
            YY(m,:) = a(m)*prod(repmat(X,numel(xj(1:m-1)),1)-repmat(xj(1:m-1),1,numel(X)),1);
        end
        Ytemp = zeros(size(X,1),size(X,2),class(x));
        for j = 1:numel(X)
            [~,I] = sort(abs(YY(:,j)));
            for m = 1:numel(a)
                Ytemp(j) = Ytemp(j) + YY(I(m),j);
            end
        end
        Y(counter,:) = Ytemp;
    else
        Y(counter,:) = a(1)*ones(size(X,1),size(X,2),class(x));
        for m = 2:numel(a)
            Y(counter,:) = Y(counter,:) + a(m)*prod(repmat(X,numel(xj(1:m-1)),1)-repmat(xj(1:m-1),1,numel(X)),1);
        end
    end
end
if plotPolynomials
    plot(X,Y)
end