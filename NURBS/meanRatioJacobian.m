function J_r = meanRatioJacobian(nurbs,xi)

d_p = nurbs.d_p;
switch d_p
    case 2
        [~,dXdxi,dXdeta] = evaluateNURBSvec(nurbs,xi,1);
        J_1 = getJacobian({dXdxi,dXdeta},NaN,d_p,false);
        J_r = 2*J_1./(norm2(dXdxi).^2 + norm2(dXdeta).^2);
    case 3
        [~,dXdxi,dXdeta,dXdzeta] = evaluateNURBSvec(nurbs,xi,1);
        J_1 = getJacobian({dXdxi,dXdeta,dXdzeta},NaN,d_p,false);
        J_r = 2*J_1./(norm2(dXdxi).^2 + norm2(dXdeta).^2 + norm2(dXdzeta).^2);
    otherwise
        error('Not implemented')
end