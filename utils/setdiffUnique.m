function Z = setdiffUnique(X,Y)

Z = zeros(size(X));
Yunique = unique(X);
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