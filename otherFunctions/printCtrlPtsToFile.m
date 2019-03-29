function printCtrlPtsToFile(coeffs, filename)

sz = size(coeffs);
M = zeros(sz(2)*sz(3),3);
for i = 1:sz(2)
    for j = 1:sz(3)
        index = i + (j-1)*sz(2);
        M(index,:) = coeffs(1:3,i,j);
    end
end
dlmwrite(filename,M,'precision',16)