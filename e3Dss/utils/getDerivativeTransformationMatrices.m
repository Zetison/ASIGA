function [J_e,J_ei,J_s,J_si,J_1,J_2] = getDerivativeTransformationMatrices(theta,phi,r)

J_e = zeros(3,3,length(theta),class(theta));
% Jacobian of spherical coordinates to cartesian for unit vectors e
J_e(1,1,:) = sin(theta).*cos(phi);
J_e(1,2,:) = sin(theta).*sin(phi);
J_e(1,3,:) = cos(theta);

J_e(2,1,:) = cos(theta).*cos(phi);
J_e(2,2,:) = cos(theta).*sin(phi);
J_e(2,3,:) = -sin(theta);

J_e(3,1,:) = -sin(phi);
J_e(3,2,:) = cos(phi);
J_e(3,3,:) = 0;


if nargout > 1
    % Inverse of J_e
    J_ei = zeros(3,3,length(theta),class(theta));
    
    J_ei(1,1,:) = sin(theta).*cos(phi);
    J_ei(1,2,:) = cos(theta).*cos(phi);
    J_ei(1,3,:) = -sin(phi);

    J_ei(2,1,:) = sin(theta).*sin(phi);
    J_ei(2,2,:) = cos(theta).*sin(phi);
    J_ei(2,3,:) = cos(phi);

    J_ei(3,1,:) = cos(theta);
    J_ei(3,2,:) = -sin(theta);
    J_ei(3,3,:) = 0;
end

if nargout > 2
    % Jacobian of spherical coordinates to cartesian coordinates
    J_s = zeros(3,3,length(theta),class(theta));
    
    J_s(1,1,:) = sin(theta).*cos(phi);
    J_s(1,2,:) = r.*cos(theta).*cos(phi);
    J_s(1,3,:) = -r.*sin(theta).*sin(phi);

    J_s(2,1,:) = sin(theta).*sin(phi);
    J_s(2,2,:) = r.*cos(theta).*sin(phi);
    J_s(2,3,:) = r.*sin(theta).*cos(phi);

    J_s(3,1,:) = cos(theta);
    J_s(3,2,:) = -r.*sin(theta);
    J_s(3,3,:) = 0;
end

if nargout > 3
    % Jacobian of spherical coordinates to cartesian coordinates
    J_si = zeros(3,3,length(theta),class(theta));
    
    J_si(1,1,:) = sin(theta).*cos(phi);
    J_si(1,2,:) = sin(theta).*sin(phi);
    J_si(1,3,:) = cos(theta);

    J_si(2,1,:) = cos(theta).*cos(phi)./r;
    J_si(2,2,:) = cos(theta).*sin(phi)./r;
    J_si(2,3,:) = -sin(theta)./r;

    J_si(3,1,:) = -sin(phi)./r; % the factor 1/sin(theta) is not taken into account here
    J_si(3,2,:) = cos(phi)./r; % the factor 1/sin(theta) is not taken into account here
    J_si(3,3,:) = 0;
end

if nargout > 4
    % Jacobian of spherical coordinates to cartesian coordinates
    J_1 = zeros(3,3,length(theta),class(theta));
    
    J_1(1,1,:) = 0;
    J_1(1,2,:) = cos(theta).*cos(phi);
    J_1(1,3,:) = -sin(phi); % the factor sin(theta) is not taken into account here

    J_1(2,1,:) = 0;
    J_1(2,2,:) = cos(theta).*sin(phi);
    J_1(2,3,:) = cos(phi); % the factor sin(theta) is not taken into account here

    J_1(3,1,:) = 0; 
    J_1(3,2,:) = -sin(theta);
    J_1(3,3,:) = 0;
    for i = 1:length(theta)
        J_1(:,:,i) = J_1(:,:,i)*J_si(:,:,i);
    end
end

if nargout > 5
    % Jacobian of spherical coordinates to cartesian coordinates
    J_2 = zeros(3,3,length(theta),class(theta));
    
    J_2(1,1,:) = 0;
    J_2(1,2,:) = -sin(theta).^2.*cos(phi); % an additional factor sin(theta) is included to compensate for the missing factor above
    J_2(1,3,:) = -cos(theta).*sin(phi);

    J_2(2,1,:) = 0;
    J_2(2,2,:) = -sin(theta).^2.*sin(phi); % an additional factor sin(theta) is included to compensate for the missing factor above
    J_2(2,3,:) = cos(theta).*cos(phi);

    J_2(3,1,:) = 0;         
    J_2(3,2,:) = -sin(theta).*cos(theta); % an additional factor sin(theta) is included to compensate for the missing factor above
    J_2(3,3,:) = 0;
    for i = 1:length(theta)
        J_2(:,:,i) = J_2(:,:,i)*J_si(:,:,i);
    end
end

