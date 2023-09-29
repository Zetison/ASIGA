function nurbs = convertToNonLinearParam(nurbs,dirs)
if nargin < 2
    dirs = 3;
else
    error('not implemented')
end
Eps = 1e4*eps;
for patch = 1:numel(nurbs)
    sz = size(nurbs{patch}.coeffs);
    n_zeta = nurbs{patch}.number(3);
    p_zeta = nurbs{patch}.degree(3);
    Zeta = nurbs{patch}.knots{3};
    uniqueKnots = unique(Zeta);
    indices = [];
    j = 1;
    for i = 1:length(uniqueKnots)
        zeta = uniqueKnots(i);
        mm = length(find(abs(Zeta - zeta) < Eps));  
        if mm >= p_zeta
            indices = [indices, j];
        end
        if i == 1
            j = j + mm - 1;
        else
            j = j + mm;
        end
    end
    pairs = [indices(1:end-1); indices(2:end)];
    for i = 1:sz(2)
        for j = 1:sz(3)
            for l = 1:size(pairs,2)
                x_start = nurbs{patch}.coeffs(1:3,i,j,pairs(1,l));
                x_end = nurbs{patch}.coeffs(1:3,i,j,pairs(2,l));
                s = linspace(0,1,pairs(2,l)-pairs(1,l)+1);
                nurbs{patch}.coeffs(1:3,i,j,pairs(1,l):pairs(2,l)) = x_start.*(1-s) + x_end.*s;
            end
        end
    end
end
    