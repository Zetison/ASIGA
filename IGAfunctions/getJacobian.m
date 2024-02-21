function [J_1,crossProd] = getJacobian(R,pts,d_p,using_R)
if nargin < 4
    using_R = true;
end

if using_R
    J1 = R{2}(:,:,1)*pts;
end
switch d_p
    case 3
        if using_R
            J2 = R{3}(:,:,1)*pts;
            J_1 = dot(J1,cross(J2,R{4}(:,:,1)*pts,2),2);
        else
            J_1 = dot(R{1},cross(R{2},R{3},2),2);
        end
        if nargout == 2
            if using_R
                crossProd = cross(J1,J2,2);
            else
                crossProd = cross(R{1},R{2},2);
            end
        end
    case 2
        if using_R
            J2 = R{3}(:,:,1)*pts;
            if size(J1,2) == 2
                crossProd = J1(:,1).*J2(:,2) - J2(:,1).*J1(:,2);
            else
                crossProd = cross(J1,J2,2);
            end
        else
            if size(R{1},2) == 2
                crossProd = R{1}(:,1).*R{2}(:,2) - R{2}(:,1).*R{1}(:,2);
            else
                crossProd = cross(R{1},R{2},2);
            end
        end
        J_1 = norm2(crossProd);
    case 1
        if using_R
            J_1 = norm2(J1); 
        else
            J_1 = norm2(R{1}); 
        end
end