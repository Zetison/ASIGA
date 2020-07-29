function dv = evaluateNURBS_deriv2(nurbs, parm_pt, derivative)
error('Depricated, use evaluateNURBS() instead')

switch nurbs.type
    case '3Dvolume'
        [~, dvdxi, dvdeta, dvdzeta] = evaluateNURBS_deriv(nurbs, parm_pt);
    case {'3Dsurface', '2Dsurface'}
        [~, dvdxi, dvdeta] = evaluateNURBS_deriv(nurbs, parm_pt);
    case {'2Dcurve', '1Dnurbs'}
        [~, dvdxi] = evaluateNURBS_deriv(nurbs, parm_pt);
end

switch derivative
    case 'xi'
        dv = dvdxi;
    case 'eta'
        dv = dvdeta;
    case 'zeta'
        dv = dvdzeta;
end