function dv = evaluateNURBS_2ndDeriv2(nurbs, parm_pt, derivative)
error('Depricated, use evaluateNURBS() instead')
switch nurbs.type
    case {'3Dsurface', '2Dsurface'}
        [~, ~, ~, d2vdxi2, d2vdxideta,d2vdeta2] = evaluateNURBS_2ndDeriv(nurbs, parm_pt);
    case {'1Dnurbs', '2Dcurve', '3Dcurve'}
        [~, ~, d2vdxi2] = evaluateNURBS_2ndDeriv(nurbs, parm_pt);
end

switch derivative
    case 'xi'
        dv = d2vdxi2;
    case 'eta'
        dv = d2vdeta2;
    case 'xieta'
        dv = d2vdxideta;
end