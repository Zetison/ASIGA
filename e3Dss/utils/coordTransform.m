function [r_m, theta_m, phi_m, A] = coordTransform(X, d_vec)
% This function first orthogonally transform the Cartesian coordinates in v
% to a cartesian coordinate system with k_vec/norm(k_vec) as the third unit
% vector. Finally, the spherical coordinates are computed in this new
% coordinate system.

[X_m, A] = orthogonalTransform(X, d_vec);

x_m = X_m(:,1);
y_m = X_m(:,2);
z_m = X_m(:,3);

%% Transform Cartesian coordinates to spherical coordinates
r_m = sqrt(x_m.^2 + y_m.^2 + z_m.^2);
theta_m = acos(z_m./r_m);
theta_m(r_m < eps) = 0;
phi_m = atan2(y_m, x_m);