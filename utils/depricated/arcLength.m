function I = arcLength(nurbs,a,b,parm_pt,dir)
error('Depricated, use NURBSarcLength instead')
 
switch dir
    case 1
        fun = @(pt) integrand(nurbs,dir,[pt, parm_pt(2)]);
    case 2
        fun = @(pt) integrand(nurbs,dir,[parm_pt(1), pt]);
    case 'eta'
        fun = @(pt) integrand(nurbs,dir,[parm_pt(1), pt, parm_pt(2)]);
end
 
n = 50;
[W1D,Q1D] = gaussianQuadNURBS(n); 
I = 0;
for i = 1:n
    I = I + fun((b-a)/2*(1+Q1D(i))+a)*W1D(i);
end
I = I*(b-a)/2;