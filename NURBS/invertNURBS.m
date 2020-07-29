function xi = invertNURBS(nurbs,x_fixed,Eps)
% Point inversion using bisection method
if isrow(x_fixed)
    x_fixed = x_fixed.';
end
if iscell(nurbs)
    warning('This routine is not written for multiple patches just yet. It should be implemented with an initial box search before using Newton''s method (using bounding boxes over each patch)')
end
d_p = nurbs.d_p;
switch d_p
    case 1
        F = @(xi) FandJacobian(xi,x_fixed,nurbs,d_p);
        xi = 0.5;
        xi = newtonsMethodND(F,NaN,xi,100,Eps,[zeros(d_p,1),ones(d_p,1)]);
    otherwise
        error('This routine should use an initial box search followed by Newton''s method')
end

function F = FandJacobian(xi,x_fixed,nurbs,d_p)

switch d_p
    case 1
        [X,dXdxi] = evaluateNURBS(nurbs,xi,1);
        F{1} = X.'-x_fixed;
        F{2} = dXdxi.';
end