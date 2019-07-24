function I = integrand(nurbs,dir,parm_pt)
 
 
switch dir
    case 1
        [~, I] = evaluateNURBS_deriv(nurbs, parm_pt);
    case {2,'eta'}
        [~, ~, I] = evaluateNURBS_deriv(nurbs, parm_pt);
end
 
I = norm(I);