function xi = invertNACA2(s,l_lm,delta_m,x_m,g,dg,dSdxi,dSdeta,s2)
xi2 = invertNACA2_2(s2,l_lm,delta_m,x_m,g,dg);
xi = zeros(size(s));
indices = s > s2;
if ~(numel(indices) == 1 && ~indices) && ~isempty(s(indices))
    xi(indices) = invertNACA2_2(s(indices),l_lm,delta_m,x_m,g,dg);
end
indices = s <= s2;
if ~(numel(indices) == 1 && ~indices) && ~isempty(s(indices))
    xi(indices) = invertNACA2_(s(indices),g,dg,dSdxi,dSdeta,s2,xi2);
end



function xi = invertNACA2_(s,g,dg,dSdxi,dSdeta,s2,xi2)
[~,integrand,arc] = NACAarcLength2(1,g,dg,dSdxi,dSdeta);
totLength = arc(xi2);
xi = zeros(size(s));
xi(1) = newtonsMethod(@(xi)arc(xi)-s(1)*totLength/s2, integrand,0.001,100,1e2*eps,[0,1]);
for i = 2:numel(s)
    [xi(i),itr] = newtonsMethod(@(xi)arc(xi)-s(i)*totLength/s2, integrand,xi(i-1)+1e-4,100,10*eps,[0,1]);
%     itr
end