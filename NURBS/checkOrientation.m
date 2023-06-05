function leftHandedOrientFound = checkOrientation(nurbs,npts)
if nargin < 2
    npts = 10;
end
xi = linspace2(0,1,npts);
[XI,ETA,ZETA] = ndgrid(xi,xi,xi);
leftHandedOrientFound = false;
for patch = 1:numel(nurbs)
    d_p = nurbs{patch}.d_p;
    if d_p == 3
        [~,dXdxi,dXdeta,dXdzeta] = evaluateNURBS(nurbs{patch},[XI(:),ETA(:),ZETA(:)],1);
        J_1 = sum(dXdxi.*cross(dXdeta,dXdzeta,2),2);
        if any(J_1 < 0)
            leftHandedOrientFound = true;
            return
        end
    else
        warning('Orientation is only controlled for d_p == 3\n')
    end
end