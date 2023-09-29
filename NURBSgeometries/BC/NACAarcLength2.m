function [s,integrand,arc] = NACAarcLength2(xi,g,dg,dSdxi,dSdeta)
dr1 = @(xi) dSdxi(xi,g(xi)) + dSdeta(xi,g(xi)).*repmat(dg(xi),1,3);

integrand = @(xi) norm2(dr1(xi.')).';
arc = @(xi) integral(integrand,0,xi);
s = arc(xi);