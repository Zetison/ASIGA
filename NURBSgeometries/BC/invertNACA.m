function xi = invertNACA(s,s2,t)

xi = zeros(size(s));
xi2 = sqrt(s2);
indices = s > s2;
if ~(numel(indices) == 1 && ~indices) && ~isempty(s(indices))
    xi(indices) = sqrt(s(indices));
end
indices = s <= s2;
if ~(numel(indices) == 1 && ~indices) && ~isempty(s(indices))
    xi(indices) = invertNACA_(s(indices),s2,xi2,t);
end


function xi = invertNACA_(s,s2,xi2,t)
df = @(xi) getNACA(xi,t,1);
integrand = @(xi) sqrt(4*xi.^2 + df(xi).^2);
arc = @(xi) NACAarcLength(xi,t);
totLength = arc(xi2);
xi = zeros(size(s));
xi(1) = newtonsMethod(@(xi)arc(xi)-s(1)*totLength/s2, integrand,0.001,100,1e2*eps,[0,1]);
for i = 2:numel(s)
    xi(i) = newtonsMethod(@(xi)arc(xi)-s(i)*totLength/s2, integrand,xi(i-1)+1e-4,100,10*eps,[0,1]);
end