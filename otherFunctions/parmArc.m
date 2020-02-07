function controlPts = parmArc(Xi,totAngle)
if Xi(end) ~= 1
    Xi = Xi/Xi(end);
end
n_xi = numel(Xi)-3;
controlPts = zeros(3,n_xi);
controlPts(:,1)   = [1, 0, 1];

counter = 1;
theta_prev = 0;

i = 4;
while i < length(Xi)
    P1 = controlPts(1:2,counter);
    
    i_prev = i;
    while sum(Xi(i) == Xi) < 2
        i = i + 1;
    end
    theta = Xi(i)*totAngle;
    P3 = [cos(theta); sin(theta)];
    t1 = (1-P1(1)*P3(1)-P1(2)*P3(2))/(P1(2)*P3(1)-P1(1)*P3(2));
    P2 = P1 + t1*[P1(2);-P1(1)];
    w2 = cos((theta-theta_prev)/2);
    
    controlPtsTemp = zeros(3,3);
    controlPtsTemp(:,1) = [P1; 1];
    controlPtsTemp(:,2) = [P2; w2];
    controlPtsTemp(:,3) = [P3; 1];
    Xi_temp = [0,0,0,1,1,1];
    nurbs_temp = createNURBSobject(controlPtsTemp,{Xi_temp});
    
    nurbs_temp = insertKnotsInNURBS(nurbs_temp,linspace2(0,1,i-i_prev));
    controlPts(:,counter+1:counter+2+i-i_prev) = nurbs_temp.coeffs(:,2:end);
    counter = counter + 2+i-i_prev;

    theta_prev = theta;
    i = i + 2;
end
