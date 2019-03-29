function [xi_rep, I] = getRepeatedKnots(Xi,p_xi)
I = [];
xi_rep = [];
grev = aveknt(Xi, p_xi+1);
for i = 1:numel(grev)
    xi = grev(i);
    if numel(find(abs(xi - Xi) < 10*eps)) >= p_xi
        I = [I, i];
        xi_rep = [xi_rep, xi];
    end
end