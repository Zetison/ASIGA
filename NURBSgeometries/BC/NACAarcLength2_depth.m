function [s,integrand,arc] = NACAarcLength2_depth(xi,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,part)

[g,dg,dSdxi,dSdeta] = getDepthRudderFunctions(b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,part(end-4:end));

dr1 = @(xi) dSdxi(xi,g(xi)) + dSdeta(xi,g(xi)).*repmat(dg(xi),3,1);

integrand = @(xi) norm2(dr1(xi).').';
arc = @(xi) integral(integrand,0,xi);
s = arc(xi);