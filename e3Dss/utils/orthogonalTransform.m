function [X_m, A] = orthogonalTransform(X, d_vec)
% This function orthogonally transform the Cartesian coordinates in v
% to a cartesian coordinate system with d_vec/norm(d_vec) as the third unit
% vector.

if isrow(d_vec)
    d_vec = d_vec.';
end
n_vec = d_vec/norm(d_vec);

%% Create orthonormal basis (e_x_m, e_y_m, n_vec)
e1 = zeros(3,1,class(X));
e1(1) = 1;

e_x_m = e1 - dot(e1,n_vec)*n_vec;
if sum(abs(e_x_m)) < 1e-15
    e_x_m = zeros(3,1,class(X));
    e_x_m(2) = 1;
end
e_x_m = e_x_m/norm(e_x_m);

e_y_m = cross(n_vec,e_x_m);

%% Create transformation matrix
A = [e_x_m, e_y_m, n_vec];

%% Do transformations
X_m = (A\X.').';