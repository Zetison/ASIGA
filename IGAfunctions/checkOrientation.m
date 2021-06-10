function checkOrientation(nurbs)

leftHandedOrientFound = false;
for patch = 1:numel(nurbs)
    d_p = nurbs{patch}.d_p;
    if d_p == 3
        [~,dXdxi,dXdeta,dXdzeta] = evaluateNURBS(nurbs{patch},0.5*ones(1,d_p),1);
        J_1 = sum(dXdxi.*cross(dXdeta,dXdzeta,2),2);
        if J_1 < 0
            fprintf('Jacobian is negative in patch %d\n', patch)
            leftHandedOrientFound = true;
        end
    else
        warning('Orientation is only controlled for d_p == 3\n')
    end
end
if ~leftHandedOrientFound
    fprintf('Great! No left handed parametrizations were found!\n')
end
