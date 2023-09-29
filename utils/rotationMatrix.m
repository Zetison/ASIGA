function R = rotationMatrix(theta, rotAxis)

R = zeros(3,3,numel(theta));
for i = 1:numel(theta)
    if isnumeric(rotAxis)
        u = rotAxis;
        u = u/norm(u);
        %https://en.wikipedia.org/wiki/Rotation_matrix
        R(:,:,i) = [cos(theta(i))+u(1)^2*(1-cos(theta(i))), u(1)*u(2)*(1-cos(theta(i)))-u(3)*sin(theta(i)), u(1)*u(3)*(1-cos(theta(i)))+u(2)*sin(theta(i));
             u(2)*u(1)*(1-cos(theta(i)))+u(3)*sin(theta(i)), cos(theta(i))+u(2)^2*(1-cos(theta(i))), u(2)*u(3)*(1-cos(theta(i)))-u(1)*sin(theta(i));
             u(3)*u(1)*(1-cos(theta(i)))-u(2)*sin(theta(i)), u(3)*u(2)*(1-cos(theta(i)))+u(1)*sin(theta(i)), cos(theta(i))+u(3)^2*(1-cos(theta(i)))];
    else
        switch rotAxis 
            case 'Xaxis'
                R(:,:,i) = [1, 0,           0;
                            0, cos(theta(i)), -sin(theta(i));
                            0, sin(theta(i)),  cos(theta(i))];
            case 'Yaxis'
                R(:,:,i) = [cos(theta(i)),  0,  sin(theta(i));
                            0,           1,  0;
                            -sin(theta(i)), 0,  cos(theta(i))];
            case 'Zaxis'
                R(:,:,i) = [cos(theta(i)),  -sin(theta(i)), 0;
                            sin(theta(i)),  cos(theta(i)),  0;
                            0,           0,           1];
        end
    end
end