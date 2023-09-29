function Z = setdiffUnique(X,Y,Eps)

if nargin < 3
    Eps = 1e-10;
end
Z = zeros(size(X));
Yunique = uniquetol(X,Eps);
counter = 1;
for i = 1:numel(Yunique)
    xi = Yunique(i);
    mX = knotRepetitions(Y,xi);
    mY = knotRepetitions(X,xi);
    if mX < mY
        mDiff = mY - mX;
        Z(counter:counter+mDiff-1) = xi;
        counter = counter + mDiff;
    end
end
Z(counter:end) = [];