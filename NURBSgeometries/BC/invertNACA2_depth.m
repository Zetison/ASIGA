function xi = invertNACA2_depth(ss,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,part,s2)
xi2 = invertNACA2_depth_2(s2,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,part);
xi = zeros(size(ss));
indices = ss > s2;
if ~(numel(indices) == 1 && ~indices) && ~isempty(ss(indices))
    xi(indices) = invertNACA2_depth_2(ss(indices),b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,part);
end
indices = ss <= s2;
if ~(numel(indices) == 1 && ~indices) && ~isempty(ss(indices))
    xi(indices) = invertNACA2_depth_(ss(indices),b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,part,s2,xi2);
end



function xi = invertNACA2_depth_(ss,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,part,s2,xi2)
[~,integrand,arc] = NACAarcLength2_depth(1,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,part);
totLength = arc(xi2);
xi = zeros(size(ss));
xi(1) = newtonsMethod(@(xi)arc(xi)-ss(1)*totLength/s2, integrand,0.001,100,1e2*eps,[0,1]);
for i = 2:numel(ss)
    xi(i) = newtonsMethod(@(xi)arc(xi)-ss(i)*totLength/s2, integrand,xi(i-1)+1e-4,100,10*eps,[0,1]);
end