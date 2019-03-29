function R = rotationMatrix(theta, rotAxis)

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