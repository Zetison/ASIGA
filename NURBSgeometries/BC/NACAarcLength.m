function s = NACAarcLength(xi,t)

df = @(xi) getNACA(xi,t,1);
integrand = @(xi) sqrt(4*xi.^2 + df(xi).^2);
arc = @(xi) integral(integrand,0,xi);
s = arc(xi);