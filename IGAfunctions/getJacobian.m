function [J_1,crossProd] = getJacobian(R,pts,d_p)

J1 = R{2}*pts;
switch d_p
    case 3
        J2 = R{3}*pts;
        J_1 = dot(J1,cross(J2,R{4}*pts,2),2);
        if nargout == 2
            crossProd = cross(J1,J2,2);
        end
    case 2
        J2 = R{3}*pts;
        crossProd = cross(J1,J2,2);
        J_1 = norm2(crossProd);
    case 1
        J_1 = norm2(J1); 
end