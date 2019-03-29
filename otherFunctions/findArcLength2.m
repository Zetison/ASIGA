function I = findArcLength2(nurbs,a,b,parm_pt,dir)
 
switch dir
    case 1
        fun = @(pt) integrand(nurbs,dir,[pt, parm_pt(2)]);
    case 2
        fun = @(pt) integrand(nurbs,dir,[parm_pt(1), pt]);
end
 
n = 50;
[W1D,Q1D] = gaussianQuadNURBS(n); 
I = 0;
for i = 1:n
    I = I + fun((b-a)/2*(1+Q1D(i))+a)*W1D(i);
end
I = I*(b-a)/2;


 
function I = integrand(nurbs,dir,parm_pt)
 
 
switch dir
    case 1
        [~, I] = evaluateNURBS_deriv(nurbs, parm_pt);
    case 2
        [~, ~, I] = evaluateNURBS_deriv(nurbs, parm_pt);
end
 
I = norm(I);