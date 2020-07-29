function vargout = evalNURBSloc(I,xi,degree,knots,weights,U,n)

R = NURBSbasis(I,xi,degree,knots,weights,n);
d_p = numel(knots);
vargout = cell(1,d_p+1);
vargout{1} = R{1}*U;
for i = 1:d_p
    for j = 1:n
        vargout{i+1}(:,:,j) = R{i+1}(:,:,j)*U;
    end
end
