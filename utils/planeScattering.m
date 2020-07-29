function data = planeScattering(X, newOptions)

% The function e3Dss (exact 3D scattering solutions) computes the solution
% to scattering problems on multilayered spherical shells impinged by a 
% plane wave or a wave due to a point source.
%
% e3Dss(X) computes  the default example evaluated at the nx3 vector X. In
% particular, the rigid scattering of a unit sphere impinged by a plane
% wave traveling along the z-axis is simulated at f=1kHz using c_f=1500m/s
% as the speed of sound.
%
% Adding the additional parameter 'options', e3Dss(X,options), enables 
% changes to the problem to be computed. In particular, options is a struct
% containing pair-wise options of the form:  
% options = struct('omega', 1e3, 'R_o', 5);
% The available parameters are the following:
% d_vec:	(default: [0;0;1]) is the direction of the plane incident wave.
% omega:	(default: 2*pi*1e3) is the angular frequency of the problem.
%          	The function accepts vector input.
% R_i:      Array of the inner radii of the spherical shells (default: [])
% R_o:      Array of the outer radii of the spherical shells (default: 1)
% P_inc:    Amplitude of incident wave at the origin. P_inc
%           can be given as a function handle P_inc(omega) where omega is the
%           angular frequency (default: 1) .
% E:        Array of the Youngs moduli of the spherical shells (default: []) 
% nu:   	Array of the Poisson ratios of the spherical shells (default: [])
% rho_s:    Array of mass densities of the spherical shells (default: []) 
% rho_f:    Array of mass densities of the fluid layers (default: []) 
% c_f:      Speed of sound in the fluid layers
%
% Author: Jon Vegard Venås
% E-mail: jon.venas@ntnu.no
% Release: 1
% Release date: 1/7/2017

options = struct('h',   	0.02,  ...      % Thickness of solid layer
                 'omega',   2*pi*1e3, ...   % Angular frequency
                 'theta1', 	0,   ...     	% Incident angle
                 'rhof1', 	1000,    ...    % Mass density of upper fluid domain
                 'rhof2', 	1.2,    ...    	% Mass density of lower fluid domain
                 'P_inc',   1,    ...     	% Amplitude of incident wave
                 'E',       210e9,   ...    % Youngs modulus for solid layers
                 'nu',      0.3,   ...    	% Poisson ratio for solid layers
                 'cf2',     340,   ...    	% Speed of sound in lower fluid domain
                 'rhos',    7850,   ...     % Mass densities of solid layer
                 'explicit',true, ...
                 'cf1',     1500);          % Speed of sound in upper fluid domain
                 
                  
if nargin > 1
    newOptionFields = fieldnames(newOptions);
    for j = 1:numel(newOptionFields)
        options.(newOptionFields{j}) = newOptions.(newOptionFields{j});
    end
end

P_inc = options.P_inc;
omega_arr = options.omega;
h = options.h;
rhos = options.rhos;
rhof1 = options.rhof1;
cf1 = options.cf1;
rhof2 = options.rhof2;
cf2 = options.cf2;
E = options.E;
nu = options.nu;

K = E/(3*(1-2*nu));
G = E/(2*(1+nu));
cs1 = sqrt((3*K+4*G)/3/rhos);
cs2 = sqrt(G/rhos);
a_arr = omega_arr/cs1;
b_arr = omega_arr/cs2;
k1_arr = omega_arr/cf1;
k2_arr = omega_arr/cf2;

if length(omega_arr) > 1
    theta1 = options.theta1;
    theta2 = asin(cs1/cf1*sin(theta1));
    gamma2 = asin(cs2/cf1*sin(theta1));
    theta3 = asin(cf2/cf1*sin(theta1));
    UU = zeros(6,numel(omega_arr));
%     for j = 1:length(omega_arr)
    parfor j = 1:length(omega_arr)
        omega = omega_arr(j);
        a = a_arr(j);
        b = b_arr(j);
        k1 = k1_arr(j);
        k2 = k2_arr(j);
        LHS = [-k1/rhof1/omega^2*cos(theta1),   a*cos(theta2),                                              -a*cos(theta2),                                             b*sin(gamma2),                                  b*sin(gamma2),                                  0; 
               1,                               a^2*(-K-G*(cos(2*theta2)+1/3)),                             a^2*(-K-G*(cos(2*theta2)+1/3)),                             -G*b^2*sin(2*gamma2),                           G*b^2*sin(2*gamma2)                             0; 
               0,                               -a^2*sin(2*theta2),                                         a^2*sin(2*theta2),                                          b^2*cos(2*gamma2)                               b^2*cos(2*gamma2)                               0; 
               0,                               -a^2*sin(2*theta2)*exp(-1i*a*h*cos(theta2)),                a^2*sin(2*theta2)*exp(1i*a*h*cos(theta2)),               	b^2*cos(2*gamma2)*exp(-1i*b*h*cos(gamma2)),    	b^2*cos(2*gamma2)*exp(1i*b*h*cos(gamma2)),      0; 
               0,                               a*cos(theta2)*exp(-1i*a*h*cos(theta2)),                     -a*cos(theta2)*exp(1i*a*h*cos(theta2)),                  	b*sin(gamma2)*exp(-1i*b*h*cos(gamma2)),        	b*sin(gamma2)*exp(1i*b*h*cos(gamma2)),          k2/rhof2/omega^2*cos(theta3)*exp(1i*k2*h*cos(theta3)); 
               0,                               a^2*(-K-G*(cos(2*theta2)+1/3))*exp(-1i*a*h*cos(theta2)),	a^2*(-K-G*(cos(2*theta2)+1/3))*exp(1i*a*h*cos(theta2)),  	-G*b^2*sin(2*gamma2)*exp(-1i*b*h*cos(gamma2)),  G*b^2*sin(2*gamma2)*exp(1i*b*h*cos(gamma2)),	exp(1i*k2*h*cos(theta3))];
        if isa(P_inc, 'function_handle')
        	RHS = [-P_inc(omega)*k1/rhof1/omega^2*cos(theta1); -P_inc(omega); 0; 0; 0; 0];
        else
        	RHS = [-P_inc*k1/rhof1/omega^2*cos(theta1); -P_inc; 0; 0; 0; 0];
        end

        Pinv = diag(1./max(abs(LHS)));
        UU(:,j) = diag(Pinv).*((LHS*Pinv)\RHS);
    end
    x = X{1}(:,1);
    z = X{1}(:,2);
    npts = size(X{1},1);
    data(1).P = repmat(UU(1,:),npts,1).*exp(1i*(x*sin(theta1)+z*cos(theta1))*k1_arr);
    data(1).dPdx = 1i*sin(theta1)*repmat(UU(1,:).*k1_arr,npts,1).*exp(1i*(x*sin(theta1)+z*cos(theta1))*k1_arr);
    data(1).dPdz = 1i*cos(theta1)*repmat(UU(1,:).*k1_arr,npts,1).*exp(1i*(x*sin(theta1)+z*cos(theta1))*k1_arr);
    x = X{2}(:,1);
    z = X{2}(:,2);
    npts = size(X{2},1);
    phi = repmat(UU(2,:),npts,1).*exp(1i*(x*sin(theta2)+z*cos(theta2))*a_arr) ...
         +repmat(UU(3,:),npts,1).*exp(1i*(x*sin(theta2)-z*cos(theta2))*a_arr);
    psi = repmat(UU(4,:),npts,1).*exp(1i*(x*sin(gamma2)+z*cos(gamma2))*b_arr) ...
         +repmat(UU(5,:),npts,1).*exp(1i*(x*sin(gamma2)-z*cos(gamma2))*b_arr);
    % u_x = d(phi)/dx - d(psi)/dz
    data(1).u_x = 1i*sin(theta2)*repmat(UU(2,:).*a_arr,npts,1).*exp(1i*(x*sin(theta2)+z*cos(theta2))*a_arr) ...
         +1i*sin(theta2)*repmat(UU(3,:).*a_arr,npts,1).*exp(1i*(x*sin(theta2)-z*cos(theta2))*a_arr) ...
         -1i*cos(gamma2)*repmat(UU(4,:).*b_arr,npts,1).*exp(1i*(x*sin(gamma2)+z*cos(gamma2))*b_arr) ...
         +1i*cos(gamma2)*repmat(UU(5,:).*b_arr,npts,1).*exp(1i*(x*sin(gamma2)-z*cos(gamma2))*b_arr);
    du_xdx =-sin(theta2)^2*repmat(UU(2,:).*a_arr.^2,npts,1).*exp(1i*(x*sin(theta2)+z*cos(theta2))*a_arr) ...
            -sin(theta2)^2*repmat(UU(3,:).*a_arr.^2,npts,1).*exp(1i*(x*sin(theta2)-z*cos(theta2))*a_arr) ...
            +0.5*sin(2*gamma2)*repmat(UU(4,:).*b_arr.^2,npts,1).*exp(1i*(x*sin(gamma2)+z*cos(gamma2))*b_arr) ...
            -0.5*sin(2*gamma2)*repmat(UU(5,:).*b_arr.^2,npts,1).*exp(1i*(x*sin(gamma2)-z*cos(gamma2))*b_arr);
    du_xdz =-0.5*sin(2*theta2)*repmat(UU(2,:).*a_arr.^2,npts,1).*exp(1i*(x*sin(theta2)+z*cos(theta2))*a_arr) ...
            +0.5*sin(2*theta2)*repmat(UU(3,:).*a_arr.^2,npts,1).*exp(1i*(x*sin(theta2)-z*cos(theta2))*a_arr) ...
            +cos(gamma2)^2*repmat(UU(4,:).*b_arr.^2,npts,1).*exp(1i*(x*sin(gamma2)+z*cos(gamma2))*b_arr) ...
            +cos(gamma2)^2*repmat(UU(5,:).*b_arr.^2,npts,1).*exp(1i*(x*sin(gamma2)-z*cos(gamma2))*b_arr);
    % u_z = d(phi)/dz + d(psi)/dx
    data(1).u_z = 1i*cos(theta2)*repmat(UU(2,:).*a_arr,npts,1).*exp(1i*(x*sin(theta2)+z*cos(theta2))*a_arr) ...
         -1i*cos(theta2)*repmat(UU(3,:).*a_arr,npts,1).*exp(1i*(x*sin(theta2)-z*cos(theta2))*a_arr) ...
         +1i*sin(gamma2)*repmat(UU(4,:).*b_arr,npts,1).*exp(1i*(x*sin(gamma2)+z*cos(gamma2))*b_arr) ...
         +1i*sin(gamma2)*repmat(UU(5,:).*b_arr,npts,1).*exp(1i*(x*sin(gamma2)-z*cos(gamma2))*b_arr);
    du_zdz =-cos(theta2)^2*repmat(UU(2,:).*a_arr.^2,npts,1).*exp(1i*(x*sin(theta2)+z*cos(theta2))*a_arr) ...
            -cos(theta2)^2*repmat(UU(3,:).*a_arr.^2,npts,1).*exp(1i*(x*sin(theta2)-z*cos(theta2))*a_arr) ...
            -0.5*sin(2*gamma2)*repmat(UU(4,:).*b_arr.^2,npts,1).*exp(1i*(x*sin(gamma2)+z*cos(gamma2))*b_arr) ...
            +0.5*sin(2*gamma2)*repmat(UU(5,:).*b_arr.^2,npts,1).*exp(1i*(x*sin(gamma2)-z*cos(gamma2))*b_arr);
    du_zdx =-0.5*sin(2*theta2)*repmat(UU(2,:).*a_arr.^2,npts,1).*exp(1i*(x*sin(theta2)+z*cos(theta2))*a_arr) ...
            +0.5*sin(2*theta2)*repmat(UU(3,:).*a_arr.^2,npts,1).*exp(1i*(x*sin(theta2)-z*cos(theta2))*a_arr) ...
            -sin(gamma2)^2*repmat(UU(4,:).*b_arr.^2,npts,1).*exp(1i*(x*sin(gamma2)+z*cos(gamma2))*b_arr) ...
            -sin(gamma2)^2*repmat(UU(5,:).*b_arr.^2,npts,1).*exp(1i*(x*sin(gamma2)-z*cos(gamma2))*b_arr);
    data(1).sigma_xx = (K+4*G/3)*du_xdx + (K-2*G/3)*du_zdz;
    data(1).sigma_yy = (K-2*G/3)*du_xdx + (K-2*G/3)*du_zdz;
    data(1).sigma_zz = (K-2*G/3)*du_xdx + (K+4*G/3)*du_zdz;
    data(1).sigma_xz = G*(du_xdz + du_zdx);
    x = X{3}(:,1);
    z = X{3}(:,2);
    npts = size(X{3},1);
    data(2).P = repmat(UU(6,:),npts,1).*exp(1i*(x*sin(theta3)-z*cos(theta3))*k2_arr);
    data(2).dPdx = 1i*sin(theta3)*repmat(UU(6,:).*k2_arr,npts,1).*exp(1i*(x*sin(theta3)-z*cos(theta3))*k2_arr);
    data(2).dPdz = -1i*cos(theta3)*repmat(UU(6,:).*k2_arr,npts,1).*exp(1i*(x*sin(theta3)-z*cos(theta3))*k2_arr);
else
    omega = omega_arr;
    theta1arr = options.theta1;
    a = a_arr;
    b = b_arr;
    k1 = k1_arr;
    k2 = k2_arr;
    P1 = zeros(size(theta1arr));
    Phi1 = zeros(size(theta1arr));
    Phi2 = zeros(size(theta1arr));
    Psi1 = zeros(size(theta1arr));
    Psi2 = zeros(size(theta1arr));
    P2 = zeros(size(theta1arr));
    theta1arr = options.theta1;
    for i = 1:length(theta1arr)
        theta1 = theta1arr(i);
        theta2 = asin(cs1/cf1*sin(theta1));
        gamma2 = asin(cs2/cf1*sin(theta1));
        theta3 = asin(cf2/cf1*sin(theta1));


        LHS = [-k1/rhof1/omega^2*cos(theta1),   a*cos(theta2),                                              -a*cos(theta2),                                             b*sin(gamma2),                                  b*sin(gamma2),                                  0; 
               1,                               a^2*(-K-G*(cos(2*theta2)+1/3)),                             a^2*(-K-G*(cos(2*theta2)+1/3)),                             -G*b^2*sin(2*gamma2),                           G*b^2*sin(2*gamma2)                             0; 
               0,                               -a^2*sin(2*theta2),                                         a^2*sin(2*theta2),                                          b^2*cos(2*gamma2)                               b^2*cos(2*gamma2)                               0; 
               0,                               -a^2*sin(2*theta2)*exp(-1i*a*h*cos(theta2)),                a^2*sin(2*theta2)*exp(1i*a*h*cos(theta2)),               	b^2*cos(2*gamma2)*exp(-1i*b*h*cos(gamma2)),    	b^2*cos(2*gamma2)*exp(1i*b*h*cos(gamma2)),      0; 
               0,                               a*cos(theta2)*exp(-1i*a*h*cos(theta2)),                     -a*cos(theta2)*exp(1i*a*h*cos(theta2)),                  	b*sin(gamma2)*exp(-1i*b*h*cos(gamma2)),        	b*sin(gamma2)*exp(1i*b*h*cos(gamma2)),          k2/rhof2/omega^2*cos(theta3)*exp(1i*k2*h*cos(theta3)); 
               0,                               a^2*(-K-G*(cos(2*theta2)+1/3))*exp(-1i*a*h*cos(theta2)),	a^2*(-K-G*(cos(2*theta2)+1/3))*exp(1i*a*h*cos(theta2)),  	-G*b^2*sin(2*gamma2)*exp(-1i*b*h*cos(gamma2)),  G*b^2*sin(2*gamma2)*exp(1i*b*h*cos(gamma2)),	exp(1i*k2*h*cos(theta3))];

        RHS = [-P_inc*k1/rhof1/omega^2*cos(theta1); -P_inc; 0; 0; 0; 0];

        Pinv = diag(1./max(abs(LHS)));
        UU = diag(Pinv).*((LHS*Pinv)\RHS);
        P1(i) = UU(1);
        Phi1(i) = UU(2);
        Phi2(i) = UU(3);
        Psi1(i) = UU(4);
        Psi2(i) = UU(5);
        P2(i) = UU(end);
    end
    data.P1 = P1;
    data.Phi1 = Phi1;
    data.Phi2 = Phi2;
    data.Psi1 = Psi1;
    data.Psi2 = Psi2;
    data.P2 = P2;
end
