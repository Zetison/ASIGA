function parm_pt = pointInversion(nurbs,x_fixed,newEpsilon)
% Point inversion using bisection method

switch nurbs.type
    case '1Dcurve'
        xi_min = 0;
        xi_max = 1;

        parm_pt = 0.5;
        x = evaluateNURBS(nurbs, parm_pt);

        while abs(x_fixed - x) > newEpsilon
            if x > x_fixed
                xi_max = parm_pt;
            else
                xi_min = parm_pt;
            end
            parm_pt = (xi_max+xi_min)/2;
            x = evaluateNURBS(nurbs, parm_pt);
        end
    case '3Dvolume'
        f = @(parm_pt) norm(x_fixed-evaluateNURBS(nurbs, parm_pt));
        options = optimset('tolx',1e-13);
        parm_pt = fminsearchbnd(f,[0.5,0.5,0.5], [0,0,0], [1,1,1], options);
%         
%         f = @(zeta) norm(x_fixed-evaluateNURBS(nurbs, [0, 1, zeta]));
% %         options = optimset();
%         parm_pt = fminsearchbnd(f,0.5, 0, 1, options);
end
% keyboard