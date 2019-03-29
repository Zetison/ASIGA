function xi = invertNACA2_2(s,l_lm,delta_m,x_m,g,dg)

f = @(xi,s) x_m - (l_lm.*xi.^2+g(xi)*delta_m.*(1-xi.^2)) - (x_m - g(0)*delta_m - s*(l_lm-g(0)*delta_m));
df = @(xi,s) -(2*l_lm.*xi+dg(xi).*delta_m.*(1-xi.^2) - 2*g(xi).*delta_m.*xi);

xi = zeros(size(s));
xi(1) = newtonsMethod(@(xi)f(xi,s(1)), @(xi)df(xi,s(1)),0.1,100,20*eps,[0,1]);
for i = 2:numel(s)
    [xi(i),itr2] = newtonsMethod(@(xi)f(xi,s(i)), @(xi)df(xi,s(i)),xi(i-1)+1e-4,100,20*eps,[0,1]);
%     itr2
end