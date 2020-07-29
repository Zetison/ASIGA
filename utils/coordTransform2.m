function [v_t, A] = coordTransform2(v, k_vec, x_0)
% This function first orthogonally transform the Cartesian coordinates in v
% to a cartesian coordinate system with k_vec/norm(k_vec) as the third unit
% vector. Finally, the the function translate the origin to x_0

if isrow(k_vec)
    k_vec = k_vec.';
end
n_vec = k_vec/norm(k_vec);

%% Create orthonormal basis (e_x_m, e_y_m, n_vec)
e1 = [1; 0; 0];

e_x_t = e1 - dot(e1,n_vec)*n_vec;
if sum(abs(e_x_t)) < 1e-15
    e_x_t = [0; 1; 0];
end
e_x_t = e_x_t/norm(e_x_t);

e_y_t = cross(n_vec,e_x_t);

%% Create transformation matrix
A = [e_x_t, e_y_t, n_vec];

%% Do transformations
v_t = (A\v.').';

x_t = v_t(:,1);
y_t = v_t(:,2);
z_t = v_t(:,3);

v_t = [x_t y_t z_t] + x_0;