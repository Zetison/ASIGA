function xi = invertNACA2_depth_2(ss,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,part)

[g,dg] = getDepthRudderFunctions(b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,part(end-4:end));

f = @(xi,ss) x_d - (l_ld.*xi.^2+g(xi).*(delta_d+(l_ud-l_ld).*xi.^2)) - (x_d - ss*l_ld);
df = @(xi) -(2*l_ld.*xi+dg(xi).*(delta_d+(l_ud-l_ld).*xi.^2)+2*g(xi).*(l_ud-l_ld).*xi);

xi = zeros(size(ss));
xi(1) = newtonsMethod(@(xi)f(xi,ss(1)), @(xi)df(xi),0.1,100,20*eps,[0,1]);
for i = 2:numel(ss)
    xi(i) = newtonsMethod(@(xi)f(xi,ss(i)), @(xi)df(xi),xi(i-1)+1e-4,100,20*eps,[0,1]);
end