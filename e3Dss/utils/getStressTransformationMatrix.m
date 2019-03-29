function D = getStressTransformationMatrix(theta,phi,type)

D = zeros(6,6,length(theta),class(theta));
switch type
    case 1
        D(1,1,:) = sin(theta).^2.*cos(phi).^2;
        D(1,2,:) = sin(theta).^2.*sin(phi).^2;
        D(1,3,:) = cos(theta).^2;
        D(1,4,:) = sin(2*theta).*sin(phi);
        D(1,5,:) = sin(2*theta).*cos(phi);
        D(1,6,:) = sin(theta).^2.*sin(2*phi);
        D(2,1,:) = cos(theta).^2.*cos(phi).^2;
        D(2,2,:) = cos(theta).^2.*sin(phi).^2;
        D(2,3,:) = sin(theta).^2;
        D(2,4,:) = -sin(2*theta).*sin(phi);
        D(2,5,:) = -sin(2*theta).*cos(phi);
        D(2,6,:) = cos(theta).^2.*sin(2*phi);
        D(3,1,:) =  sin(phi).^2;
        D(3,2,:) = cos(phi).^2;
        
        D(3,6,:) = -sin(2*phi);
        D(4,1,:) = -0.5*cos(theta).*sin(2*phi);
        D(4,2,:) = 0.5*cos(theta).*sin(2*phi);
        
        D(4,4,:) = -sin(theta).*cos(phi);
        D(4,5,:) = sin(theta).*sin(phi);
        D(4,6,:) = cos(theta).*cos(2*phi);
        D(5,1,:) = -0.5*sin(theta).*sin(2*phi);
        D(5,2,:) = 0.5*sin(theta).*sin(2*phi);
        
        D(5,4,:) = cos(theta).*cos(phi);
        D(5,5,:) = -cos(theta).*sin(phi);
        D(5,6,:) = sin(theta).*cos(2*phi);
        D(6,1,:) = 0.5*sin(2*theta).*cos(phi).^2;
        D(6,2,:) = 0.5*sin(2*theta).*sin(phi).^2;
        D(6,3,:) = -0.5*sin(2*theta);
        D(6,4,:) = cos(2*theta).*sin(phi);
        D(6,5,:) = cos(2*theta).*cos(phi);
        D(6,6,:) = 0.5*sin(2*theta).*sin(2*phi);
          
    case 2 % Compute the inverse matrix

        D(1,1,:) = sin(theta).^2.*cos(phi).^2;
        D(1,2,:) = cos(theta).^2.*cos(phi).^2;
        D(1,3,:) = sin(phi).^2;
        D(1,4,:) = -cos(theta).*sin(2*phi);
        D(1,5,:) = -sin(theta).*sin(2*phi);
        D(1,6,:) = sin(2*theta).*cos(phi).^2;
        D(2,1,:) = sin(theta).^2.*sin(phi).^2;
        D(2,2,:) = cos(theta).^2.*sin(phi).^2;
        D(2,3,:) = cos(phi).^2;
        D(2,4,:) = cos(theta).*sin(2*phi);
        D(2,5,:) = sin(theta).*sin(2*phi);
        D(2,6,:) = sin(2*theta).*sin(phi).^2;
        D(3,1,:) = cos(theta).^2;
        D(3,2,:) = sin(theta).^2;

        D(3,6,:) = -sin(2*theta);
        D(4,1,:) = 0.5*sin(2*theta).*sin(phi);
        D(4,2,:) = -0.5*sin(2*theta).*sin(phi);

        D(4,4,:) = -sin(theta).*cos(phi);
        D(4,5,:) = cos(theta).*cos(phi);
        D(4,6,:) = sin(phi).*cos(2*theta);
        D(5,1,:) = 0.5*sin(2*theta).*cos(phi);
        D(5,2,:) = -0.5*sin(2*theta).*cos(phi);

        D(5,4,:) = sin(theta).*sin(phi);
        D(5,5,:) = -cos(theta).*sin(phi);
        D(5,6,:) = cos(2*theta).*cos(phi);
        D(6,1,:) = 0.5*sin(theta).^2.*sin(2*phi);
        D(6,2,:) = 0.5*cos(theta).^2.*sin(2*phi);
        D(6,3,:) = -0.5*sin(2*phi);
        D(6,4,:) = cos(theta).*cos(2*phi);
        D(6,5,:) = cos(2*phi).*sin(theta);
        D(6,6,:) = 0.5*sin(2*theta).*sin(2*phi);
end
