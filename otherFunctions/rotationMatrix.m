function R = rotationMatrix(theta, rotAxis)

if isnumeric(rotAxis)
    u = rotAxis;
    u = u/norm(u);
    %https://en.wikipedia.org/wiki/Rotation_matrix
    R = [cos(theta)+u(1)^2*(1-cos(theta)), u(1)*u(2)*(1-cos(theta))-u(3)*sin(theta), u(1)*u(3)*(1-cos(theta))+u(2)*sin(theta);
         u(2)*u(1)*(1-cos(theta))+u(3)*sin(theta), cos(theta)+u(2)^2*(1-cos(theta)), u(2)*u(3)*(1-cos(theta))-u(1)*sin(theta);
         u(3)*u(1)*(1-cos(theta))-u(2)*sin(theta), u(3)*u(2)*(1-cos(theta))+u(1)*sin(theta), cos(theta)+u(3)^2*(1-cos(theta))];
else
    switch rotAxis 
        case 'Xaxis'
            R = [1, 0,           0;
                 0, cos(theta), -sin(theta);
                 0, sin(theta),  cos(theta)];
        case 'Yaxis'
            R = [cos(theta),  0,  sin(theta);
                 0,           1,  0;
                 -sin(theta), 0,  cos(theta)];
        case 'Zaxis'
            R = [cos(theta),  -sin(theta), 0;
                 sin(theta),  cos(theta),  0;
                 0,           0,           1];
    end
end