function dirs = findDirsNSD2(js,cps,dgdx,dgdy,dim,zeta)
if any(isinf(cps(2,:)))
    dirs = NaN;
    return
end

if nargin == 3
    dg = @(xi,r) dgdx{r}(xi);
elseif dim == 1
    dg = @(xi,r) dgdx{r}(xi,zeta);
else
    dg = @(xi,r) dgdy{r}(zeta,xi);
end

dirs = zeros(size(js));
for i = 1:size(cps,2)
    xi = cps(1,i); 
    r = cps(2,i);
    sgndgdx = sign(dg(xi,r+1));
    if i == 1
        dirs(1) = exp((4*js(1)+sgndgdx)*pi*1i/2/(r+1));
    elseif i == size(cps,2)
        dirs(end) = exp((4*js(end)+sgndgdx)*pi*1i/2/(r+1));
    else
        dirs(2*i-2:2*i-1) = exp((4*js(2*i-2:2*i-1)+sgndgdx)*pi*1i/2/(r+1));
    end
end