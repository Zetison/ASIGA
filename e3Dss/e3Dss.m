function data = e3Dss(X, newOptions)

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
% Eps:      Small parameter for series truncation. The summation are
%           terminated whenever the relative contribution of the given term is
%           less then Eps. If vector input is given for either X or omega, the
%           maximal relative contribution of the given term is compared
%           with Eps
% N_max:    Upper limit for the number of terms in the series
% prec:     Precision of the calculations (default: 'double'). Both 'sym'
%           and 'mp' are supported, with arbitrary precision altered by Digits and
%           mp.Digits respectively
% calc_farField: Set true if the far field is to be calculated for outer
%           fluid domain
% usePlaneWave: Set true for plane incident waves
% usePointChargeWave: Set true for point charge incident waves (the
%           location of the point source is given by d_vec*r_s
% r_s:      Radius to source location for point charge incident waves
% calc_p:   Calculate the scattered pressure
% calc_dpdx: Calculate the derivative of the scattered pressure fields w.r.t. the cartesian coordinate x
% calc_dpdy: Calculate the derivative of the scattered pressure fields w.r.t. the cartesian coordinate y
% calc_dpdz: Calculate the derivative of the scattered pressure fields w.r.t. the cartesian coordinate z
% calc_u_x: Calculate the x component of the displacement field
% calc_u_y: Calculate the y component of the displacement field 
% calc_u_z: Calculate the z component of the displacement field
% calc_u_r: Calculate the r component of the displacement field
% calc_u_t: Calculate the theta component of the displacement field
% calc_du_xdx: Calculate the derivative of the x component of the displacement field w.r.t. x
% calc_du_xdy: Calculate the derivative of the x component of the displacement field w.r.t. y
% calc_du_xdz: Calculate the derivative of the x component of the displacement field w.r.t. z
% calc_du_ydx: Calculate the derivative of the y component of the displacement field w.r.t. x
% calc_du_ydy: Calculate the derivative of the y component of the displacement field w.r.t. y
% calc_du_ydz: Calculate the derivative of the y component of the displacement field w.r.t. z
% calc_du_zdx: Calculate the derivative of the z component of the displacement field w.r.t. x
% calc_du_zdy: Calculate the derivative of the z component of the displacement field w.r.t. y
% calc_du_zdz: Calculate the derivative of the z component of the displacement field w.r.t. z
% calc_sigma_xx: Calculate the derivative of the xx component of the stress field (cartesian coordinates)
% calc_sigma_yy: Calculate the derivative of the yy component of the stress field (cartesian coordinates)
% calc_sigma_zz: Calculate the derivative of the zz component of the stress field (cartesian coordinates)
% calc_sigma_yz: Calculate the derivative of the yz component of the stress field (cartesian coordinates)
% calc_sigma_xz: Calculate the derivative of the xz component of the stress field (cartesian coordinates)
% calc_sigma_xy: Calculate the derivative of the xy component of the stress field (cartesian coordinates)
% calc_p_laplace: Calculate the Laplace operator of the scattered pressure fields
% calc_errorsDisplacementCondition: Calculate the errors for the displacement conditions
% calc_errorsPressureCondition: Calculate the errors for the pressure conditions
% calc_errorsHelmholtz: Calculate errors for the Helmholtz equation
% calc_errorsNavier:    Calculate errors for the Navier equation
%
% Note that only the minimal set of parameters should be given. For
% example, if numel(E) < numel(R_o) it is assumed that the inner domain is
% an elastic sphere. Moreover, to exclude the inner (fluid) domain by
% implementing sound soft boundary conditions on the inner surface, simply
% reduce numel(rho_f) and numel(c_f) correspondingly. Sound hard boundary
% conditions may be implemented in a similar fashion on the inner sphere by
% reducing numel(E), numel(nu) and numel(rho_s) correspondingly.
%
% If the solution is needed in domains other then the outer (fluid) domain,
% then X should be passed as a cell array, where the points in X{1}
% corresponds to the outer domain X{2} the first solid domain etc.
%
% The function handles the case in which both X and omega are vectors.
% 
% The following improvement will be implemented in the future:
%   - Implement a scaling in the linear system of equations such that
%   evaluations of the spherical bessel functions j_n(x) and y_n(x) are
%   replaced by j_n(x)/x^n and y_n(x)*x^(n+1), respectively. This will in
%   turn solve the problem of evaluations extending the range of the given
%   precision (10^290 is implemented for double), and thus enable
%   evaluations at even higher frequencies.
%   - Implement fluid-fluid layers, and solid-solid layers
%   - Implement polynomial expansions for thin solid layers to avoid bad
%   conditioning of matrix
%   - Implement Fenders strategy of avoiding complex matrices.
%   - Implement Fenders strategy of reducing the size of the matrix.
%   - Generalize to general input sources on both the outer and inner
%   surface.
%
% Author: Jon Vegard Venås
% E-mail: jon.venas@ntnu.no
% Release: 1
% Release date: 1/7/2017

%% Dirichlet and Robin conditions has not been implemented.
options = struct('d_vec',   [0;0;1],  ... 	% Direction of the incident wave
                 'omega',   2*pi*1e3, ...   % Angular frequency
                 'R_i',     [],   ...     	% Inner radii
                 'R_o',     1,    ...    	% Outer radii
                 'P_inc',   1,    ...     	% Amplitude of incident wave
                 'E',       [],   ...      	% Youngs modulus for solid layers
                 'nu',      [],   ...    	% Poisson ratio for solid layers
                 'rho_s',   [],   ...    	% Mass densities of solid layers
                 'rho_f',   [],   ...       % Mass densities of fluid layers
                 'c_f',     1500, ...       % Speed of sound in fluid layers
                 'Eps',     eps,  ...       % Small parameter for series truncation
                 'N_max',   inf,  ...       % Upper limit for the number of terms in the series
                 'prec',    'double', ...   % Precision of the calculations
                 'calc_farField',       false, ...  % Set true if the far field is to be calculated for outer fluid domain
                 'usePlaneWave',        true,  ...  % Set true for plane incident waves
                 'usePointChargeWave',  false, ...  % Set true for point charge incident waves
                 'r_s',                 NaN,   ...	% Radius to source location for point charge incident waves
                 'calc_p',              true,  ...  % Calculate the scattered pressure
                 'calc_dpdx',           false, ...  % Calculate the derivative of the scattered pressure fields w.r.t. the cartesian coordinate x
                 'calc_dpdy',           false, ...  % Calculate the derivative of the scattered pressure fields w.r.t. the cartesian coordinate y
                 'calc_dpdz',           false, ...  % Calculate the derivative of the scattered pressure fields w.r.t. the cartesian coordinate z
                 'calc_u_x',            false, ...  % Calculate the x component of the displacement
                 'calc_u_y',            false, ...  % Calculate the y component of the displacement
                 'calc_u_z',            false, ...  % Calculate the z component of the displacement
                 'calc_du_xdx',         false, ...  % Calculate the derivative of the x component of the displacement field w.r.t. x
                 'calc_du_xdy',         false, ...  % Calculate the derivative of the x component of the displacement field w.r.t. y
                 'calc_du_xdz',         false, ...  % Calculate the derivative of the x component of the displacement field w.r.t. z
                 'calc_du_ydx',         false, ...  % Calculate the derivative of the y component of the displacement field w.r.t. x
                 'calc_du_ydy',         false, ...  % Calculate the derivative of the y component of the displacement field w.r.t. y
                 'calc_du_ydz',         false, ...  % Calculate the derivative of the y component of the displacement field w.r.t. z
                 'calc_du_zdx',         false, ...  % Calculate the derivative of the z component of the displacement field w.r.t. x
                 'calc_du_zdy',         false, ...  % Calculate the derivative of the z component of the displacement field w.r.t. y
                 'calc_du_zdz',         false, ...  % Calculate the derivative of the z component of the displacement field w.r.t. z
                 'calc_sigma_xx',       false, ...  % Calculate the derivative of the xx component of the stress field (cartesian coordinates)
                 'calc_sigma_yy',       false, ...  % Calculate the derivative of the yy component of the stress field (cartesian coordinates)
                 'calc_sigma_zz',       false, ...  % Calculate the derivative of the zz component of the stress field (cartesian coordinates)
                 'calc_sigma_yz',       false, ...  % Calculate the derivative of the yz component of the stress field (cartesian coordinates)
                 'calc_sigma_xz',       false, ...  % Calculate the derivative of the xz component of the stress field (cartesian coordinates)
                 'calc_sigma_xy',       false, ...  % Calculate the derivative of the xy component of the stress field (cartesian coordinates)
                 'calc_p_laplace',      false, ...  % Calculate the Laplace operator of the scattered pressure fields
                 'calc_errorsDisplacementCondition', false, ... % Calculate the errors for the displacement conditions
                 'calc_errorsPressureCondition', 	 false, ... % Calculate the errors for the pressure conditions
                 'calc_errorsHelmholtz',             false, ... % Calculate errors for the Helmholtz equation
                 'calc_errorsNavier',                false);	% Calculate errors for the Navier equation
                  
if nargin > 1
    newOptionFields = fieldnames(newOptions);
    for j = 1:numel(newOptionFields)
        options.(newOptionFields{j}) = newOptions.(newOptionFields{j});
    end
end

if ~iscell(X)
    X = {X};
end
prec = options.prec;

omega = options.omega;
E = options.E;
nu = options.nu;
rho_s = options.rho_s;
if numel(E) ~= numel(nu) || numel(E) ~= numel(rho_s) || numel(nu) ~= numel(rho_s)
    error('The number of parameters in E, nu and rho_s must all be the same')
end
rho_f = options.rho_f;
R_o = options.R_o;
R_i = options.R_i;
d_vec = options.d_vec;

calc_errorsPressureCondition = options.calc_errorsPressureCondition;
calc_errorsDisplacementCondition = options.calc_errorsDisplacementCondition;
calc_errorsHelmholtz = options.calc_errorsHelmholtz;
calc_errorsNavier = options.calc_errorsNavier;
calc_p = options.calc_p || calc_errorsPressureCondition;
calc_dpdx = options.calc_dpdx || calc_errorsDisplacementCondition;
calc_dpdy = options.calc_dpdy || calc_errorsDisplacementCondition;
calc_dpdz = options.calc_dpdz || calc_errorsDisplacementCondition;
calc_u_x = options.calc_u_x;
calc_u_y = options.calc_u_y;
calc_u_z = options.calc_u_z;

calc_du_xdx = options.calc_du_xdx;
calc_du_xdy = options.calc_du_xdy;
calc_du_xdz = options.calc_du_xdz;
calc_du_ydx = options.calc_du_ydx;
calc_du_ydy = options.calc_du_ydy;
calc_du_ydz = options.calc_du_ydz;
calc_du_zdx = options.calc_du_xdx;
calc_du_zdy = options.calc_du_zdy;
calc_du_zdz = options.calc_du_zdz;
calc_sigma_xx = options.calc_sigma_xx;
calc_sigma_yy = options.calc_sigma_yy;
calc_sigma_zz = options.calc_sigma_zz;
calc_sigma_yz = options.calc_sigma_yz;
calc_sigma_xz = options.calc_sigma_xz;
calc_sigma_xy = options.calc_sigma_xy;
calc_p_laplace = options.calc_p_laplace;
calc_errors = calc_errorsDisplacementCondition || calc_errorsPressureCondition || calc_errorsHelmholtz || calc_errorsNavier;

options.calc_dpdr = calc_dpdx || calc_dpdy || calc_dpdz || calc_p_laplace;
options.calc_dpdt = options.calc_dpdr;
options.calc_d2pdr2 = calc_p_laplace || calc_errorsHelmholtz;
options.calc_d2pdt2 = calc_p_laplace || calc_errorsHelmholtz;
calcCartesianDispDerivatives = calc_du_xdx || calc_du_xdy || calc_du_xdz || calc_du_ydx || calc_du_ydy || calc_du_ydz || calc_du_zdx || calc_du_zdy || calc_du_zdz;
options.calc_u_r = calc_errorsDisplacementCondition || calcCartesianDispDerivatives || calc_u_x || calc_u_y || calc_u_z || calc_errorsNavier;
options.calc_u_t = calc_errorsDisplacementCondition || calcCartesianDispDerivatives || calc_u_x || calc_u_y || calc_u_z || calc_errorsNavier;
options.calc_du_rdr = calcCartesianDispDerivatives;
options.calc_du_rdt = calcCartesianDispDerivatives;
options.calc_du_tdr = calcCartesianDispDerivatives;
options.calc_du_tdt = calcCartesianDispDerivatives;

calcStresses =  calc_sigma_xx || calc_sigma_yy || calc_sigma_zz || calc_sigma_yz || calc_sigma_xz ...
                             || calc_sigma_xy || calc_errorsNavier || calc_errorsPressureCondition;
options.calc_sigma_rr = calcStresses;
options.calc_sigma_tt = calcStresses;
options.calc_sigma_pp = calcStresses;
options.calc_sigma_rt = calcStresses;
options.calc_errorsNavier = calc_errorsNavier;

if isrow(omega)
    options.omega = omega.';
end
if options.usePointChargeWave 
    if isnan(options.r_s)
        options.r_s = 2*R_o(1);
    elseif abs(options.r_s) < R_o(1)
        error('r_s must satisfy |r_s| > R_o(1)')
    end
    warning(['It is not implemented an efficient routine to evaluate Equation (D.7). ' ...
             'The built in MATLAB routine is used for this purpose.'])
end
if options.usePointChargeWave
    options.usePlaneWave = false;
end
if ~(options.usePointChargeWave || options.usePlaneWave)
    error('Only plane waves and point charge waves are implemented for the incident wave.')
end
if numel(E) < numel(R_o)
    SHBC = true;
    ESBC = false;
    SSBC = false;  
else
    SHBC = false;
    if length(R_o) > length(R_i) % Inner domain is a solid sphere
        ESBC = true;
        SSBC = false;  
    else
        ESBC = false;
        if length(R_o) == length(rho_f)
            SSBC = true;
        else
            SSBC = false;        
        end
    end
end
options.ESBC = ESBC;
options.SHBC = SHBC;
options.SSBC = SSBC;

if omega(1) == 0 && options.usePointChargeWave
    error('This case has no unique solution')
end
    

%% Coordinate transformation
% Due to symmetry, we can do a coordinate transformation such that we only
% need to compute the solution for the special case k_vec = k*[0, 0, 1].
r = cell(size(X));
theta = cell(size(X));
phi = cell(size(X));
for j = 1:length(X)
    if ~isempty(X{j})
        [r{j}, theta{j}, phi{j}, A] = coordTransform(X{j}, d_vec);
        r{j} = r{j}.';
        theta{j} = theta{j}.';
        phi{j} = phi{j}.';    
    end
end

%% Compute the solution with d_vec = [0, 0, 1]
data_0 = e3Dss_0(r, theta, options);

%% Coordinate transformation (Transform back to original Cartesian coordinate system)

% Allocate memory
nFreqs = length(omega);

M = length(options.R_o);
if SHBC || ESBC || SSBC
    data(M).p = [];
else
    data(M+1).p = [];
end
m = 1;
for j = 1:length(X)
    if ~isempty(r{j})
        if mod(j,2)
            if calc_p
                data(m).p = zeros(size(X{j},1),nFreqs,prec);
            end
            if calc_dpdx
                data(m).dpdx = zeros(size(X{j},1),nFreqs,prec);
            end
            if calc_dpdy
                data(m).dpdy = zeros(size(X{j},1),nFreqs,prec);
            end
            if calc_dpdz
                data(m).dpdz = zeros(size(X{j},1),nFreqs,prec);
            end
            if calc_p_laplace
                data(m).p_laplace = zeros(size(X{j},1),nFreqs,prec);
            end
        else
            if calc_u_x
                data(m).u_x = zeros(size(X{j},1),nFreqs,prec);
            end
            if calc_u_y
                data(m).u_y = zeros(size(X{j},1),nFreqs,prec);
            end
            if calc_u_z
                data(m).u_z = zeros(size(X{j},1),nFreqs,prec);
            end
        end
    end
    if ~mod(j,2)
        m = m + 1;
    end
end

m = 1;
for j = 1:length(X)
    if calc_errors
        if mod(j,2)
            data(m).v_fluid = X{j};
        else
            data(m).v_solid = X{j};
        end
    end
    if ~isempty(r{j})
        Theta = repmat(theta{j},nFreqs,1);
        Phi = repmat(phi{j},nFreqs,1);
        R = repmat(r{j},nFreqs,1);
        indices = logical(R < eps);
        if mod(j,2)
            if calc_dpdx || calc_dpdy || calc_dpdz || calc_p_laplace
                dpdr = data_0(m).dpdr;
                dpdt = sin(Theta).*data_0(m).dpdt; % rescale dpdt
                dpdX_m = cell(3,1);

                dpdX_m{1} = dpdr.*sin(Theta).*cos(Phi) + dpdt.*cos(Theta).*cos(Phi)./R;
                dpdX_m{1}(indices) = 0;
                
                dpdX_m{2} = dpdr.*sin(Theta).*sin(Phi) + dpdt.*cos(Theta).*sin(Phi)./R;
                dpdX_m{2}(indices) = 0;

                dpdX_m{3} = dpdr.*cos(Theta) - dpdt.*sin(Theta)./R;
                dpdX_m{3}(indices) = dpdr(indices);

                dpdX = cell(3,1);
                dpdX{1} = zeros(size(X{j},1),nFreqs);
                dpdX{2} = zeros(size(X{j},1),nFreqs);
                dpdX{3} = zeros(size(X{j},1),nFreqs);
                for ii = 1:3
                    for jj = 1:3
                        dpdX{ii} = dpdX{ii} + A(ii,jj)*dpdX_m{jj}.';
                    end
                end
                if calc_dpdx
                    data(m).dpdx = dpdX{1};
                end
                if calc_dpdy
                    data(m).dpdy = dpdX{2};
                end
                if calc_dpdz
                    data(m).dpdz = dpdX{3};
                end
                dpdt = data_0(m).dpdt; % use scaled version of dpdt for the laplace operator
                if calc_p_laplace
                    d2pdr2 = data_0(m).d2pdr2;
                    d2pdt2 = data_0(m).d2pdt2;

                    temp = 2./R.*dpdr + 1./R.^2.*cos(Theta).*dpdt + 1./R.^2.*d2pdt2;                    
                    temp(indices) = 0;
                    data(m).p_laplace = d2pdr2 + temp;
                    data(m).p_laplace = data(m).p_laplace.';
                end
            end
            if calc_p
                data(m).p = data_0(m).p.';
            end
        else
            if calcStresses
                % Transform the stresses in the spherical coordinate system to the 
                % Cartesian coordinate system
                sigma_m = cell(6,1);
                sigma_m{1} = data_0(m).sigma_rr;
                sigma_m{2} = data_0(m).sigma_tt;
                sigma_m{3} = data_0(m).sigma_pp;
                sigma_m{4} = zeros(nFreqs,size(X{j},1),prec);
                sigma_m{5} = zeros(nFreqs,size(X{j},1),prec);
                sigma_m{6} = data_0(m).sigma_rt;
                
                D = getStressTransformationMatrix(theta{j},phi{j},2);
                sigma_X_m = cell(6,1);
                for ii = 1:6
                    sigma_X_m{ii} = zeros(nFreqs,size(X{j},1),prec);
                end
                sigma_X = sigma_X_m;
                for ii = 1:6
                    for jj = 1:6
                        D_kl = D(ii, jj, :);
                        D_kl = D_kl(:);
                        sigma_X_m{ii} = sigma_X_m{ii} + repmat(D_kl.',nFreqs,1).*sigma_m{jj};
                    end
                end
                if ESBC && m == M
                    sigma_X_m{1}(indices) = sigma_m{1}(indices);
                    sigma_X_m{2}(indices) = sigma_m{2}(indices);
                    sigma_X_m{3}(indices) = sigma_m{3}(indices);
                    sigma_X_m{4}(indices) = 0;
                    sigma_X_m{5}(indices) = 0;
                    sigma_X_m{6}(indices) = 0;
                end

                alpha = zeros(3,prec);
                I = eye(3);
                for ii = 1:3
                    for jj = 1:3
                        alpha(ii,jj) = dot(I(:,ii), A(:,jj));
                    end
                end

                vgt = [1 1;
                       2 2;
                       3 3;
                       2 3;
                       1 3;
                       1 2];
                vgtinv = [1 6 5;
                          6 2 4;
                          5 4 3];
                for vgtIdx = 1:6
                    for ii = 1:3
                        for jj = 1:3
                            sigma_X{vgtIdx} = sigma_X{vgtIdx} + alpha(vgt(vgtIdx,1),ii)*alpha(vgt(vgtIdx,2),jj)*sigma_X_m{vgtinv(ii,jj)};
                        end
                    end
                end
                for ii = 1:6
                    sigma_X{ii} = sigma_X{ii}.';
                end
                if calc_sigma_xx
                    data(m).sigma_xx = sigma_X{1};
                end
                if calc_sigma_yy
                    data(m).sigma_yy = sigma_X{2};
                end
                if calc_sigma_zz
                    data(m).sigma_zz = sigma_X{3};
                end
                if calc_sigma_yz
                    data(m).sigma_yz = sigma_X{4};
                end
                if calc_sigma_xz
                    data(m).sigma_xz = sigma_X{5};
                end
                if calc_sigma_xy
                    data(m).sigma_xy = sigma_X{6};
                end
            end
            if calcCartesianDispDerivatives
                % Transform the derivatives in the spherical coordinate system to the 
                % Cartesian coordinate system
                u_r = data_0(m).u_r;
                u_t = data_0(m).u_t;
                
                du_m = cell(3,3);
                du_m{1,1} = data_0(m).du_rdr;
                du_m{1,2} = data_0(m).du_rdt.*sin(Theta); % rescale du_rdt
                du_m{1,3} = zeros(nFreqs,size(X{j},1),prec);
                du_m{2,1} = data_0(m).du_tdr.*sin(Theta); % rescale du_tdr
                du_m{2,2} = data_0(m).du_tdt;
                du_m{2,3} = zeros(nFreqs,size(X{j},1),prec);
                du_m{3,1} = zeros(nFreqs,size(X{j},1),prec);
                du_m{3,2} = zeros(nFreqs,size(X{j},1),prec);
                du_m{3,3} = zeros(nFreqs,size(X{j},1),prec);
                
                [J_e,~,~,J_si,J_1,J_2] = getDerivativeTransformationMatrices(theta{j},phi{j},r{j});
                du_X_m = cell(3,3);
                for ii = 1:3
                    for jj = 1:3
                        du_X_m{ii,jj} = zeros(nFreqs,size(X{j},1),prec);
                    end
                end
                du_X = du_X_m;
                temp = du_X_m;
                for ii = 1:3
                    for jj = 1:3
                        for ll = 1:3
                            J_si_lj = J_si(ll, jj, :);
                            J_si_lj = J_si_lj(:);
                            temp{ii,jj} = temp{ii,jj} + du_m{ii,ll}.*repmat(J_si_lj.',nFreqs,1);
                        end
                    end
                end
                for ii = 1:3
                    for jj = 1:3
                        for ll = 1:3
                            J_e_lj = J_e(ll, ii, :); % transpose of J_e ...
                            J_e_lj = J_e_lj(:);
                            du_X_m{ii,jj} = du_X_m{ii,jj} + repmat(J_e_lj.',nFreqs,1).*temp{ll,jj};
                        end
                    end
                end
                for ii = 1:3
                    for jj = 1:3
                        J_1_ij = J_1(ii, jj, :);
                        J_1_ij = J_1_ij(:);
                        du_X_m{ii,jj} = du_X_m{ii,jj} + repmat(J_1_ij.',nFreqs,1).*u_r;
                        
                        J_2_ij = J_2(ii, jj, :);
                        J_2_ij = J_2_ij(:);
                        du_X_m{ii,jj} = du_X_m{ii,jj} + repmat(J_2_ij.',nFreqs,1).*u_t;
                    end
                end
                if ESBC && m == M
                    du_X_m{1,1}(indices) = data_0(m).du_rdr(indices);
                    du_X_m{2,2}(indices) = data_0(m).du_rdt(indices);
                    du_X_m{3,3}(indices) = data_0(m).du_tdr(indices);
                    for ii = 1:3
                        for jj = 1:3
                            if ii ~= jj
                                du_X_m{ii,jj}(indices) = 0;
                            end
                        end
                    end
                end

                for ii = 1:3
                    for jj = 1:3
                        temp{ii,jj}(:) = 0;
                    end
                end
                Ainv = inv(A);
                for ii = 1:3
                    for jj = 1:3
                        for ll = 1:3
                            temp{ii,jj} = temp{ii,jj} + A(ii, ll)*du_X_m{ll,jj};
                        end
                    end
                end
                for ii = 1:3
                    for jj = 1:3
                        for ll = 1:3
                            du_X{ii,jj} = du_X{ii,jj} + temp{ii, ll}*Ainv(ll,jj);
                        end
                    end
                end
                for ii = 1:3
                    for jj = 1:3
                        du_X{ii,jj} = du_X{ii,jj}.';
                    end
                end
                
                if calc_du_xdx
                    data(m).du_xdx = du_X{1,1};
                end
                if calc_du_xdy
                    data(m).du_xdy = du_X{1,2};
                end
                if calc_du_xdz
                    data(m).du_xdz = du_X{1,3};
                end
                if calc_du_ydx
                    data(m).du_ydx = du_X{2,1};
                end
                if calc_du_ydy
                    data(m).du_ydy = du_X{2,2};
                end
                if calc_du_ydz
                    data(m).du_ydz = du_X{2,3};
                end
                if calc_du_zdx
                    data(m).du_zdx = du_X{3,1};
                end
                if calc_du_zdy
                    data(m).du_zdy = du_X{3,2};
                end
                if calc_du_zdz
                    data(m).du_zdz = du_X{3,3};
                end
            end
            if calc_u_x || calc_u_y || calc_u_z
                u_X_m = cell(3,1);
                u_r = data_0(m).u_r;
                u_t = data_0(m).u_t.*sin(Theta); % rescale u_t
                u_X_m{1} = u_r.*sin(Theta).*cos(Phi) + u_t.*cos(Theta).*cos(Phi);
                u_X_m{2} = u_r.*sin(Theta).*sin(Phi) + u_t.*cos(Theta).*sin(Phi);
                u_X_m{3} = u_r.*cos(Theta) - u_t.*sin(Theta);
                if ESBC && m == M
                    u_X_m{1}(indices) = 0;
                    u_X_m{2}(indices) = 0;
                    u_X_m{3}(indices) = u_r(indices);
                end
                u_X = cell(3,1);
                u_X{1} = zeros(size(X{j},1),nFreqs,prec);
                u_X{2} = zeros(size(X{j},1),nFreqs,prec);
                u_X{3} = zeros(size(X{j},1),nFreqs,prec);
                
                for ii = 1:3
                    for jj = 1:3
                        u_X{ii} = u_X{ii} + A(ii,jj)*u_X_m{jj}.';
                    end
                end
                if calc_u_x
                    data(m).u_x = u_X{1};
                end
                if calc_u_y
                    data(m).u_y = u_X{2};
                end
                if calc_u_z
                    data(m).u_z = u_X{3};
                end
            end
            if calc_errorsNavier
                data(m).navier1 = data_0(m).navier1.';
                data(m).navier2 = data_0(m).navier2.';
            end
        end
    end
    if ~mod(j,2)
        m = m + 1;
    end
end
if calc_errors
    for m = 1:M+1
        if m ~= M+1
            for investigate = {'innerSurface','outerSurface'}
                if m ~= M+1 && ~((SHBC || ESBC) && m == M && strcmp(investigate{1},'innerSurface'))
                    iS = strcmp(investigate{1},'innerSurface');
                    Eps = options.Eps;
                    if SSBC && m == M && strcmp(investigate{1},'innerSurface')
                        v_f = X{2*m};
                    else
                        v_f = X{2*(m+iS)-1};
                        c_f = options.c_f(m+iS);
                    end
                    if m == M && SHBC
                        indices_f = abs(norm2(X{2*m-1}) - R_o(m)) < 10*Eps;
                    elseif SSBC && m == M && strcmp(investigate{1},'innerSurface')
                        indices_f = abs(norm2(X{2*m}) - R_i(m)) < 10*Eps;    
                        indices_s = indices_f;
                    else
                        [indices_f, indices_s] = findMatchingPoints(v_f,X{2*m},Eps);
                    end
                    v_f = v_f(indices_f,:);
                    P_inc = options.P_inc;


                    if calc_errorsDisplacementCondition && ~(SSBC && m == M && strcmp(investigate{1},'innerSurface'))
                        rho_f = options.rho_f(m+iS);
                        
                        n_x = v_f(:,1)./norm2(v_f);
                        n_y = v_f(:,2)./norm2(v_f);
                        n_z = v_f(:,3)./norm2(v_f);
                        n_x = repmat(n_x,1,length(omega));
                        n_y = repmat(n_y,1,length(omega));
                        n_z = repmat(n_z,1,length(omega));
                        
                        dpdx = data(m+iS).dpdx(indices_f,:);
                        dpdy = data(m+iS).dpdy(indices_f,:);
                        dpdz = data(m+iS).dpdz(indices_f,:);
                        if m == 1 && ~iS
                            k = omega/c_f(1);
                            k_vec = d_vec*k;
                            p_inc = P_inc*exp(1i*v_f*k_vec);
                            dpdx = dpdx + 1i*p_inc.*repmat(k_vec(1,:),size(p_inc,1),1);
                            dpdy = dpdy + 1i*p_inc.*repmat(k_vec(2,:),size(p_inc,1),1);
                            dpdz = dpdz + 1i*p_inc.*repmat(k_vec(3,:),size(p_inc,1),1);
                        end
                        if SHBC && m == M && strcmp(investigate{1},'outerSurface')
                            err_dc = max(abs(dpdx.*n_x + dpdy.*n_y + dpdz.*n_z),[],1)/P_inc;
                        else
                            u_x = data(m).u_x(indices_s,:);
                            u_y = data(m).u_y(indices_s,:);
                            u_z = data(m).u_z(indices_s,:);
                            Omega = repmat(omega,size(u_x,1),1);
                            err_dc = max(abs(     (dpdx-rho_f*Omega.^2.*u_x).*n_x ...
                                                + (dpdy-rho_f*Omega.^2.*u_y).*n_y ...
                                                + (dpdz-rho_f*Omega.^2.*u_z).*n_z),[],1)./max(abs(dpdx.*n_x + dpdy.*n_y + dpdz.*n_z),[],1);
                        end
                        if ~isfield(data(m),'err_dc') || isempty(data(m).err_dc)
                            data(m).err_dc = err_dc;
                        else
                            data(m).err_dc = max([data(m).err_dc; err_dc],[],1);
                        end
                    end

                    if calc_errorsPressureCondition
                        if ~(SHBC && m == M)
                            if ~(SSBC && m == M && strcmp(investigate{1},'innerSurface'))
                                p_tot = data(m+iS).p(indices_f,:);
                                if m == 1 && ~iS
                                    k = omega/c_f(1);
                                    k_vec = d_vec*k;
                                    p_inc = P_inc*exp(1i*v_f*k_vec);
                                    p_tot = p_tot + p_inc;
                                end
                            end
                            sigma_cart = cell(6,1);
                            sigma_cart{1} = data(m).sigma_xx(indices_s,:);
                            sigma_cart{2} = data(m).sigma_yy(indices_s,:);
                            sigma_cart{3} = data(m).sigma_zz(indices_s,:);
                            sigma_cart{4} = data(m).sigma_yz(indices_s,:);
                            sigma_cart{5} = data(m).sigma_xz(indices_s,:);
                            sigma_cart{6} = data(m).sigma_xy(indices_s,:);

                            sigma_rr = zeros(size(sigma_cart{1}),prec);
                            
                            phi_temp = atan2(v_f(:,2),v_f(:,1));
                            r_temp = sqrt(v_f(:,1).^2+v_f(:,2).^2+v_f(:,3).^2);
                            theta_temp = acos(v_f(:,3)./r_temp);
                            D = getStressTransformationMatrix(theta_temp,phi_temp,1);
                            for l = 1:6
                                D_kl = D(1, l, :);
                                D_kl = repmat(D_kl(:),1,length(omega));
                                sigma_rr = sigma_rr + D_kl.*sigma_cart{l};
                            end
                            if SSBC && m == M && strcmp(investigate{1},'innerSurface')
                                err_pc = max(abs(sigma_rr),[],1)/P_inc;
                            else
                                err_pc = max(abs(p_tot+sigma_rr),[],1)./max(abs(p_tot),[],1);
                            end
                            if ~isfield(data(m), 'err_pc') || isempty(data(m).err_pc)
                                data(m).err_pc = err_pc;
                            else
                                data(m).err_pc = max([data(m).err_pc; err_pc],[],1);
                            end
                        end
                    end
                end
            end
            if calc_errorsNavier && ~(SHBC && m == M)
                Theta = repmat(theta{2*m},nFreqs,1);
                rho_s = options.rho_s(m);
                u_r = data_0(m).u_r;
                u_t = data_0(m).u_t.*sin(Theta); % rescale u_t
                u_r = u_r.';
                u_t = u_t.';
                Omega = repmat(omega,size(u_r,1),1);
                data(m).err_navier1 = max(abs(data(m).navier1 + rho_s*Omega.^2.*u_r),[],1)./max(abs(rho_s*Omega.^2.*u_r),[],1);
                data(m).err_navier2 = max(abs(data(m).navier2 + rho_s*Omega.^2.*u_t),[],1)./max(abs(rho_s*Omega.^2.*u_t),[],1);
            end
        end
        if calc_errorsHelmholtz && ~(m == M+1 && (SSBC || SHBC || ESBC))
            p_laplace = data(m).p_laplace;
            k = omega/options.c_f(m);
            K = repmat(k,size(data(m).p,1),1);

            data(m).err_helmholtz = max(abs(p_laplace+K.^2.*data(m).p),[],1)./max(abs(K.^2.*data(m).p),[],1);
        end
    end
end
data(1).flag = data_0(1).flag;
data(1).N_eps = data_0(1).N_eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = e3Dss_0(r, theta, options)
% This function computes exact 3D scattering solutions when the axis of
% symmetry is the z-axis
%
%
% Note that dpdt, u_t, du_tdr, and du_rdt are scaled by csc(theta)

omega = options.omega;
c_f = options.c_f;
E = options.E;
R_o = options.R_o;

prec = options.prec;

SSBC = options.SSBC;
ESBC = options.ESBC;
SHBC = options.SHBC;

Eps = options.Eps;
nFreqs = length(omega);
M = length(R_o);
N_max = options.N_max;

k = omega*(1./c_f);

% Compute derived quantities
if ~(SHBC && M == 1)
    nu = options.nu;
    rho_s = options.rho_s;
    K = E./(3*(1-2*nu));
    G = E./(2*(1+nu));
    options.G = G;

    c_s_1 = sqrt((3*K+4*G)./(3*rho_s)); % longitudinal wave velocity 
    c_s_2 = sqrt(G./rho_s); % shear wave velocity
    
    a = omega*(1./c_s_1);
    b = omega*(1./c_s_2);
end


if any(omega == 0)
    computeForStaticCase = true;
else
    computeForStaticCase = false;
end
% Allocate memory
P = cell(length(r),1);
dP = cell(length(r),1);
d2P = cell(length(r),1);
if ~(SHBC && M == 1)
    Z(M).xi = cell(2,2); % Spherical Bessel functions evaluated at a*r
    Z(M).eta = cell(2,2); % Spherical Bessel functions evaluated at b*r
end
if SHBC || ESBC || SSBC
    data(M).p = [];
    Z(M).zeta = cell(2,2); % Spherical Bessel functions evaluated at k*r
else
    data(M+1).p = [];
    Z(M+1).zeta = cell(2,2); % Spherical Bessel functions evaluated at k*r
end

m = 1;
for j = 1:length(r)
    if ~isempty(r{j})
        P{j} = zeros(2,length(theta{j}),prec); 
        dP{j} = zeros(2,length(theta{j}),prec); 
        d2P{j} = zeros(2,length(theta{j}),prec);
        if mod(j,2)
            Z(m).zeta{1,1} = zeros(nFreqs-computeForStaticCase,length(theta{j}),prec);
            Z(m).zeta{1,2} = zeros(nFreqs-computeForStaticCase,length(theta{j}),prec);
            Z(m).zeta{2,1} = zeros(nFreqs-computeForStaticCase,length(theta{j}),prec);
            Z(m).zeta{2,2} = zeros(nFreqs-computeForStaticCase,length(theta{j}),prec);
            if options.calc_p
                data(m).p = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_dpdr
                data(m).dpdr = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_dpdt
                data(m).dpdt = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_d2pdr2
                data(m).d2pdr2 = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_d2pdt2
                data(m).d2pdt2 = zeros(nFreqs,length(r{j}),prec);
            end
        else
            Z(m).xi{1,1} = zeros(nFreqs-computeForStaticCase,length(theta{j}),prec);
            Z(m).xi{1,2} = zeros(nFreqs-computeForStaticCase,length(theta{j}),prec);
            Z(m).xi{2,1} = zeros(nFreqs-computeForStaticCase,length(theta{j}),prec);
            Z(m).xi{2,2} = zeros(nFreqs-computeForStaticCase,length(theta{j}),prec);

            Z(m).eta{1,1} = zeros(nFreqs-computeForStaticCase,length(theta{j}),prec);
            Z(m).eta{1,2} = zeros(nFreqs-computeForStaticCase,length(theta{j}),prec);
            Z(m).eta{2,1} = zeros(nFreqs-computeForStaticCase,length(theta{j}),prec);
            Z(m).eta{2,2} = zeros(nFreqs-computeForStaticCase,length(theta{j}),prec);
            if options.calc_u_r
                data(m).u_r = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_u_t
                data(m).u_t = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_du_rdr
                data(m).du_rdr = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_du_rdt
                data(m).du_rdt = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_du_tdr
                data(m).du_tdr = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_du_tdt
                data(m).du_tdt = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_sigma_rr
                data(m).sigma_rr = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_sigma_tt
                data(m).sigma_tt = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_sigma_pp
                data(m).sigma_pp = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_sigma_rt
                data(m).sigma_rt = zeros(nFreqs,length(r{j}),prec);
            end
            if options.calc_errorsNavier
                data(m).navier1 = zeros(nFreqs,length(r{j}),prec);
                data(m).navier2 = zeros(nFreqs,length(r{j}),prec);
            end
        end
    end
    if ~mod(j,2)
        m = m + 1;
    end
end

indices = (1:nFreqs).';
Zindices = indices;
m = 1;
if computeForStaticCase
    if isa(options.P_inc,'function_handle')
        P_inc = options.P_inc(0);
    else
        P_inc = options.P_inc;
    end
    staticIdx = find(omega == 0);
    for j = 1:length(r)
        if ~isempty(r{j})
            if mod(j,2)
                if m > 1
                    data(m).p(staticIdx,:) = P_inc*ones(1,length(r{j}),prec);
                end
            else
                A = -P_inc/(3*K);
                if options.calc_u_r
                    data(m).u_r(staticIdx,:) = A*repmat(r{j},nFreqs,1);
                end
                if options.calc_sigma_rr
                    data(m).sigma_rr(staticIdx,:) = A*ones(1,length(r{j}),prec);
                end
                if options.calc_sigma_tt
                    data(m).sigma_tt(staticIdx,:) = A*ones(1,length(r{j}),prec);
                end
                if options.calc_sigma_pp
                    data(m).sigma_pp(staticIdx,:) = A*ones(1,length(r{j}),prec);
                end
            end
        end
        if ~mod(j,2)
            m = m + 1;
        end
    end
    indices(staticIdx) = [];
    Zindices = Zindices(1:end-1);
end

% In order to avoid premeture termination of the series summation, the
% vector wasNonZero{m} contains the information about the magnitude of the
% previous contribution to the sum. In particular, if the previous
% nExtraTerms terms have relative magnitude less than Eps compared to the
% total sum, the summation will terminate (for the particular omega)
nExtraTerms = 2;
hasCnvrgd = zeros(nFreqs,nExtraTerms); % matrix of element that "has converged"
if computeForStaticCase
    hasCnvrgd(staticIdx,:) = 1;
end
if strcmp(prec,'sym')
    tiny = vpa('1e-1000'); % To avoid dividing by zero.
else
    tiny = realmin(prec); % To avoid dividing by zero.
end
countUpwards = 1;
if countUpwards
    n = zeros(1,prec);
else
    n = round(3*max(r{1})*max(k(:,1)));
end
data(1).N_eps = NaN(nFreqs,1);
flag = zeros(size(hasCnvrgd,1),1); % Program terminated successfully unless error occurs (for each frequency)
while n >= 0 && n <= N_max
    try % and hope that no spherical Bessel functions are evaluated to be too large
        omega_temp = omega(indices);
        k_temp = k(indices,:);
        if SHBC && M == 1
            a_temp = [];
            b_temp = [];
        else
            a_temp = a(indices,:);
            b_temp = b(indices,:);
        end
        
        CC = getCoeffs(n, omega_temp, a_temp, b_temp, k_temp, options);

        m = 1;
        hasCnvrgdTmp = zeros(length(indices),length(r)); % temporary hasCnvrgd matrix
        for j = 1:length(r)    
            hasCnvrgdTmp2 = ones(size(indices)); % temporary hasCnvrgd vector    
            if ~isempty(r{j})
                if countUpwards
                    [P{j}, dP{j}, d2P{j}] = legendreDerivs(n, cos(theta{j}), P{j}, dP{j}, d2P{j});
                else
                    P{j}(2,:) = legendre(n, cos(theta{j}));
                    dP{j}(2,:) = legendreDeriv(n, cos(theta{j}));
                    d2P{j}(2,:) = legendreDeriv2(n, cos(theta{j}));
                end
                if mod(j,2)
                    if m == 1
                        C = CC(:,1);
                    elseif m == M + 1
                        C = CC(:,end);
                    else
                        C = CC(:,6*(m-1):6*(m-1)+1);
                    end
                    zeta = k_temp(:,m)*r{j};
                    if ~options.calc_farField
                        if n == 0
                            Z(m).zeta{1,2} = bessel_s(n,zeta,1);
                            if m < M+1
                                Z(m).zeta{2,2} = bessel_s(n,zeta,2);
                            end
                        end
                        Z(m).zeta{1,1} = Z(m).zeta{1,2}(Zindices,:);
                        Z(m).zeta{1,2} = bessel_s(n+1,zeta,1);
                        if m < M+1
                            Z(m).zeta{2,1} = Z(m).zeta{2,2}(Zindices,:);
                            Z(m).zeta{2,2} = bessel_s(n+1,zeta,2);    
                        end
                    end
                    fluid = p_(m,n,zeta,theta{j},M,C,k_temp(:,m),P{j}(2,:), dP{j}(2,:), d2P{j}(2,:), Z(m).zeta, options);
                    if options.calc_p
                        data(m).p(indices,:) = data(m).p(indices,:) + fluid.p;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(fluid.p)./(abs(data(m).p(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_dpdr
                        data(m).dpdr(indices,:) = data(m).dpdr(indices,:) + fluid.dpdr;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(fluid.dpdr)./(abs(data(m).dpdr(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_dpdt
                        data(m).dpdt(indices,:) = data(m).dpdt(indices,:) + fluid.dpdt;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(fluid.dpdt)./(abs(data(m).dpdt(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_d2pdr2
                        data(m).d2pdr2(indices,:) = data(m).d2pdr2(indices,:) + fluid.d2pdr2;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(fluid.d2pdr2)./(abs(data(m).d2pdr2(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_d2pdt2
                        data(m).d2pdt2(indices,:) = data(m).d2pdt2(indices,:) + fluid.d2pdt2;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(fluid.d2pdt2)./(abs(data(m).d2pdt2(indices,:))+tiny) < Eps,2);
                    end
                else
                    C_indices = 6*(m-1)+2:6*m-1;
                    
                    xi = a_temp(:,m)*r{j};
                    eta = b_temp(:,m)*r{j};
                    if n == 0
                        Z(m).xi{1,2} = bessel_s(n,xi,1);
                        Z(m).eta{1,2} = bessel_s(n,eta,1);
                        if ~(ESBC && m == M)
                            Z(m).xi{2,2} = bessel_s(n,xi,2);
                            Z(m).eta{2,2} = bessel_s(n,eta,2);
                        end
                    end     
                    Z(m).xi{1,1} = Z(m).xi{1,2}(Zindices,:);
                    Z(m).eta{1,1} = Z(m).eta{1,2}(Zindices,:);
                    Z(m).xi{1,2} = bessel_s(n+1,xi,1);
                    Z(m).eta{1,2} = bessel_s(n+1,eta,1);
                    if ~(ESBC && m == M)
                        Z(m).xi{2,1} = Z(m).xi{2,2}(Zindices,:);
                        Z(m).eta{2,1} = Z(m).eta{2,2}(Zindices,:);
                        Z(m).xi{2,2} = bessel_s(n+1,xi,2); 
                        Z(m).eta{2,2} = bessel_s(n+1,eta,2);
                    end
                    if ESBC && m == M
                        A = CC(:,C_indices(1));
                        B = CC(:,C_indices(3));
                    else
                        A = CC(:,C_indices(1:2));
                        B = CC(:,C_indices(3:4));
                    end
                    options.indices = indices;
                    solid = u_(m,n,r{j},theta{j},M,A,B,xi,eta,G(m),K(m),a_temp(:,m),b_temp(:,m),...
                                P{j}(2,:),dP{j}(2,:),d2P{j}(2,:), Z(m).xi, Z(m).eta, ESBC, options);

                    if options.calc_u_r
                        data(m).u_r(indices,:) = data(m).u_r(indices,:) + solid.u_r;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.u_r)./(abs(data(m).u_r(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_u_t
                        data(m).u_t(indices,:) = data(m).u_t(indices,:) + solid.u_t;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.u_t)./(abs(data(m).u_t(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_du_rdr
                        data(m).du_rdr(indices,:) = data(m).du_rdr(indices,:) + solid.du_rdr;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.du_rdr)./(abs(data(m).du_rdr(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_du_rdt
                        data(m).du_rdt(indices,:) = data(m).du_rdt(indices,:) + solid.du_rdt;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.du_rdt)./(abs(data(m).du_rdt(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_du_tdr
                        data(m).du_tdr(indices,:) = data(m).du_tdr(indices,:) + solid.du_tdr;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.du_tdr)./(abs(data(m).du_tdr(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_du_tdt
                        data(m).du_tdt(indices,:) = data(m).du_tdt(indices,:) + solid.du_tdt;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.du_tdt)./(abs(data(m).du_tdt(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_sigma_rr
                        data(m).sigma_rr(indices,:) = data(m).sigma_rr(indices,:) + solid.sigma_rr;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.sigma_rr)./(abs(data(m).sigma_rr(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_sigma_tt
                        data(m).sigma_tt(indices,:) = data(m).sigma_tt(indices,:) + solid.sigma_tt;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.sigma_tt)./(abs(data(m).sigma_tt(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_sigma_pp
                        data(m).sigma_pp(indices,:) = data(m).sigma_pp(indices,:) + solid.sigma_pp;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.sigma_pp)./(abs(data(m).sigma_pp(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_sigma_rt
                        data(m).sigma_rt(indices,:) = data(m).sigma_rt(indices,:) + solid.sigma_rt;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.sigma_rt)./(abs(data(m).sigma_rt(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_errorsNavier
                        data(m).navier1(indices,:) = data(m).navier1(indices,:) + solid.navier1;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.navier1)./(abs(data(m).navier1(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_errorsNavier
                        data(m).navier2(indices,:) = data(m).navier2(indices,:) + solid.navier2;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.navier2)./(abs(data(m).navier2(indices,:))+tiny) < Eps,2);
                    end
                end
            end
            hasCnvrgdTmp(logical(hasCnvrgdTmp2),j) = 1;
            
            if ~mod(j,2)
                m = m + 1;
            end
        end
        if countUpwards 
            hasCnvrgd(indices,:) = [hasCnvrgd(indices,2:end), prod(hasCnvrgdTmp,2)];
            indicesPrev = indices;
            indices = find(~prod(hasCnvrgd,2));
            [~,Zindices] = ismember(indices,indicesPrev);
            if length(indices) < length(indicesPrev)
                data(1).N_eps(setdiff(indicesPrev,indices)) = n;
            end
            if isempty(indices) % every element has converged
                break;
            end
            n = n + 1;
        else
            n = n - 1;
        end
    catch        
        flag = -~prod(hasCnvrgd,2);
        warning(['The summation ended prematurely at n = ' num2str(n) ...
                 ' because a Bessel function evaluation was too large or that the' ...
                 ' global matrix was singular to working precision.'])
        break
    end
end
data(1).flag = flag;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fluid = p_(m,n,zeta,theta,M,C,k,P,dP,d2P,Z,options)
% Note that in the case of m = M+1 and zeta = 0: 
% --- dpdz =: dpdr and dpdx = dpdy = dpdt = 0
% --- nabla p =: d2pdr2, d2pdt2 := 0
% Also note that dpdt is scaled by csc(theta)

Q0 = P;
if options.calc_farField
    h_n   = repmat(1i^(-n-1)./k, 1, size(zeta,2));
    if options.calc_dpdr
        dh_n  = repmat(1i^(-n), size(zeta,1),size(zeta,2));
    end
    if options.calc_d2pdr2
        d2h_n = repmat(1i^(-n+1)*k, 1, size(zeta,2));
    end
else
    j_n = Z{1,1};
    if options.calc_dpdr
        dj_n = dbessel_s(n,zeta,1,Z);
    end
    if options.calc_d2pdr2
        d2j_n = d2bessel_s(n,zeta,1,Z);
    end
    if m == 1
        y_n = Z{2,1};
        
        h_n = j_n + 1i*y_n;
        if options.calc_dpdr
            dy_n = dbessel_s(n,zeta,2,Z);
            dh_n = dj_n + 1i*dy_n;
        end
        if options.calc_d2pdr2
            d2y_n = d2bessel_s(n,zeta,2,Z);
            d2h_n = d2j_n + 1i*d2y_n;
        end
    elseif m < M+1
        y_n = Z{2,1};
        if options.calc_dpdr
            dy_n = dbessel_s(n,zeta,2,Z);
        end
        if options.calc_d2pdr2
            d2y_n = d2bessel_s(n,zeta,2,Z);
        end
    end
end

if options.calc_dpdt
    Q1 = Q_(1,theta,P,dP,d2P,true);
end
if options.calc_d2pdt2
    Q2 = Q_(2,theta,P,dP,d2P);
end

if m == 1
    if options.calc_p
        fluid.p = C*Q0.*h_n;
    end
    if options.calc_dpdr
        fluid.dpdr = C.*k*Q0.*dh_n;
    end
    if options.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        fluid.dpdt = C*Q1.*h_n;
    end
    if options.calc_d2pdr2
        fluid.d2pdr2 = C.*k.^2*Q0.*d2h_n;
    end
    if options.calc_d2pdt2
        fluid.d2pdt2 = C*Q2.*h_n;
    end
elseif m == M+1
    if options.calc_p
        fluid.p = C*Q0.*j_n;
    end
    if options.calc_dpdr
        fluid.dpdr = C.*k*Q0.*dj_n;
        indices = logical(zeta(1,:) < eps);
        if n == 1
            fluid.dpdr(:,indices) = repmat(k/3.*C,1,sum(indices));
        else
            fluid.dpdr(:,indices) = 0;
        end
    end
    if options.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        fluid.dpdt = C*Q1.*j_n;
        indices = logical(zeta(1,:) < eps);
        fluid.dpdt(:,indices) = 0;
    end
    if options.calc_d2pdr2
        fluid.d2pdr2 = C.*k.^2*Q0.*d2j_n;
        indices = logical(zeta(1,:) < eps);
        if n == 0
            fluid.d2pdr2(:,indices) = repmat(-k.^2.*C,1,sum(indices));
        else
            fluid.d2pdr2(:,indices) = 0;
        end
    end
    if options.calc_d2pdt2
        fluid.d2pdt2 = C*Q2.*j_n;
        indices = logical(zeta(1,:) < eps);
        fluid.d2pdt2(:,indices) = 0;
    end
else
    if options.calc_p
        fluid.p = C(:,1)*Q0.*j_n + C(:,2)*Q0.*y_n;
    end
    if options.calc_dpdr
        fluid.dpdr = C(:,1).*k*Q0.*dj_n + C(:,2).*k*Q0.*dy_n;
    end
    if options.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        fluid.dpdt = C(:,1)*Q1.*j_n + C(:,2)*Q1.*y_n;
    end
    if options.calc_d2pdr2
        fluid.d2pdr2 = C(:,1).*k.^2*Q0.*d2j_n + C(:,2).*k.^2*Q0.*d2y_n;
    end
    if options.calc_d2pdt2
        fluid.d2pdt2 = C(:,1)*Q2.*j_n + C(:,2)*Q2.*y_n;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solid = u_(m,n,r,theta,M,A,B,xi,eta,G,K,a,b,P,dP,d2P,Zxi,Zeta,ESBC,options)
% Note that in the case of "ESBC" and r = 0: 
% --- u_z =: u_r and u_x = u_y = u_t = 0
% --- du_xdx =: du_rdr, du_ydy =: du_rdt, du_zdz =: du_tdr, 0 =: du_tdt
% --- sigma_11 =: sigma_rr, sigma_22 =: sigma_tt, sigma_33 =: sigma_pp, 0 =: sigma_rt
% Also note that u_t, du_tdr and du_rdt are scaled by csc(theta)

Q0 = Q_(0,theta,P,dP,d2P);
Q1s = Q_(1,theta,P,dP,d2P,true);
Q1 = sin(theta).*Q1s;
Q2 = Q_(2,theta,P,dP,d2P);

r2 = r.^2;

if options.calc_u_r
    Q0r = Q0./r;
    u_r =   A(:,1)*Q0r.*S_(1,1,n,xi,eta,Zxi) ...
          + B(:,1)*Q0r.*T_(1,1,n,eta,Zeta);
    if ESBC && m == M
        indices = logical(r < eps);
        if n == 1
            u_r(:,indices) = repmat((a.*A(:,1) - 2*b.*B(:,1))/3,1,sum(indices));
        else
            u_r(:,indices) = 0;
        end
    else
        u_r = u_r + A(:,2)*Q0r.*S_(1,2,n,xi,eta,Zxi) ...
                  + B(:,2)*Q0r.*T_(1,2,n,eta,Zeta);
    end
    solid.u_r = u_r;
end

if options.calc_u_t
    Q1sr = Q1s./r;
    u_t =   A(:,1)*Q1sr.*S_(2,1,n,xi,eta,Zxi) ...
          + B(:,1)*Q1sr.*T_(2,1,n,eta,Zeta);
    if ESBC && m == M
        u_t(:,logical(r < eps)) = 0;
    else
        u_t = u_t + A(:,2)*Q1sr.*S_(2,2,n,xi,eta,Zxi) ...
                  + B(:,2)*Q1sr.*T_(2,2,n,eta,Zeta);
    end
    solid.u_t = u_t;
end

if options.calc_du_rdr
    Q0r2 = Q0./r2;
    du_rdr =   A(:,1)*Q0r2.*S_(3,1,n,xi,eta,Zxi) ...
             + B(:,1)*Q0r2.*T_(3,1,n,eta,Zeta);
    if ESBC && m == M
        indices = logical(r < eps);
        if n == 0
            du_rdr(:,indices) = repmat(G/K*(4*a.^2-3*b.^2).*A(:,1)/9,1,sum(indices));
        elseif n == 2
            du_rdr(:,indices) = repmat(-(a.^2.*A(:,1)-3*b.^2.*B(:,1))/15,1,sum(indices));
        else
            du_rdr(:,indices) = 0;
        end
    else
        du_rdr = du_rdr + A(:,2)*Q0r2.*S_(3,2,n,xi,eta,Zxi) ...
                        + B(:,2)*Q0r2.*T_(3,2,n,eta,Zeta);
    end
    solid.du_rdr = du_rdr;
end

if options.calc_du_rdt
    Q1sr = Q1s./r;
    du_rdt =   A(:,1)*Q1sr.*S_(1,1,n,xi,eta,Zxi) ...
             + B(:,1)*Q1sr.*T_(1,1,n,eta,Zeta);
    if ESBC && m == M
        indices = logical(r < eps);
        if n == 0
            du_rdt(:,indices) = repmat(G/K*(4*a.^2-3*b.^2).*A(:,1)/9,1,sum(indices));
        elseif n == 2
            du_rdt(:,indices) = repmat(-(a.^2.*A(:,1)-3*b.^2.*B(:,1))/15,1,sum(indices));
        else
            du_rdt(:,indices) = 0;
        end
    else
        du_rdt = du_rdt + A(:,2)*Q1sr.*S_(1,2,n,xi,eta,Zxi) ...
                        + B(:,2)*Q1sr.*T_(1,2,n,eta,Zeta);
    end
    solid.du_rdt = du_rdt;
end

if options.calc_du_tdr
    Q1sr2 = Q1s./r2;
    du_tdr =   A(:,1)*Q1sr2.*S_(4,1,n,xi,eta,Zxi) ...
             + B(:,1)*Q1sr2.*T_(4,1,n,eta,Zeta);
    if ESBC && m == M
        indices = logical(r < eps);
        if n == 0
            du_tdr(:,indices) = repmat(G/K*(4*a.^2-3*b.^2).*A(:,1)/9,1,sum(indices));
        elseif n == 2
            du_tdr(:,indices) = repmat(2*(a.^2.*A(:,1)-3*b.^2.*B(:,1))/15,1,sum(indices));
        else
            du_tdr(:,indices) = 0;
        end
    else
        du_tdr = du_tdr + A(:,2)*Q1sr2.*S_(4,2,n,xi,eta,Zxi) ...
                        + B(:,2)*Q1sr2.*T_(4,2,n,eta,Zeta);
    end
    solid.du_tdr = du_tdr;
end

if options.calc_du_tdt
    Q2r1 = Q2./r;
    du_tdt =   A(:,1)*Q2r1.*S_(2,1,n,xi,eta,Zxi) ...
             + B(:,1)*Q2r1.*T_(2,1,n,eta,Zeta);
    if ESBC && m == M
        du_tdt(:,logical(r < eps)) = 0;
    else
        du_tdt = du_tdt + A(:,2)*Q2r1.*S_(2,2,n,xi,eta,Zxi) ...
                        + B(:,2)*Q2r1.*T_(2,2,n,eta,Zeta);
    end
    solid.du_tdt = du_tdt;
end

if options.calc_sigma_rr
    Q0r2 = Q0./r2;
    sigma_rr =   A(:,1)*Q0r2.*S_(5,1,n,xi,eta,Zxi) ...
               + B(:,1)*Q0r2.*T_(5,1,n,eta,Zeta);
    if ESBC && m == M
        indices = logical(r < eps);
        if n == 0
            sigma_rr(:,indices) = repmat(5*(4*a.^2-3*b.^2).*A(:,1)/30,1,sum(indices));
        elseif n == 2
            sigma_rr(:,indices) = repmat((-2*a.^2.*A(:,1) + 6*b.^2.*B(:,1))/30,1,sum(indices));
        else
            sigma_rr(:,indices) = 0;
        end
    else
        sigma_rr = sigma_rr + A(:,2)*Q0r2.*S_(5,2,n,xi,eta,Zxi) ...
                            + B(:,2)*Q0r2.*T_(5,2,n,eta,Zeta);
    end
	sigma_rr = 2*G*sigma_rr;
    
	solid.sigma_rr = sigma_rr;
end
if options.calc_sigma_tt
    Q0r2 = Q0./r2;
    Q2r2 = Q2./r2;
    sigma_tt =   A(:,1)*Q0r2.*S_(6,1,n,xi,eta,Zxi)  ...
               + B(:,1)*Q0r2.*T_(6,1,n,eta,Zeta) ...
               + A(:,1)*Q2r2.*S_(2,1,n,xi,eta,Zxi) ...
               + B(:,1)*Q2r2.*T_(2,1,n,eta,Zeta);
    if ESBC && m == M
        indices = logical(r < eps);
        if n == 0
            sigma_tt(:,indices) = repmat(5*(4*a.^2-3*b.^2).*A(:,1)/30,1,sum(indices));
        elseif n == 2
            sigma_tt(:,indices) = repmat((-2*a.^2.*A(:,1) + 6*b.^2.*B(:,1))/30,1,sum(indices));
        else
            sigma_tt(:,indices) = 0;
        end
    else
        sigma_tt = sigma_tt + A(:,2)*Q0r2.*S_(6,2,n,xi,eta,Zxi) ...
                            + B(:,2)*Q0r2.*T_(6,2,n,eta,Zeta) ...
                            + A(:,2)*Q2r2.*S_(2,2,n,xi,eta,Zxi) ...
                            + B(:,2)*Q2r2.*T_(2,2,n,eta,Zeta);
    end
	sigma_tt = 2*G*sigma_tt;
    
	solid.sigma_tt = sigma_tt;
end
if options.calc_sigma_pp
    Q0r2 = Q0./r2;
    cott_Q1r2 = cos(theta).*Q1s./r2;
    sigma_pp =   A(:,1)*Q0r2.*S_(6,1,n,xi,eta,Zxi) ...
               + B(:,1)*Q0r2.*T_(6,1,n,eta,Zeta) ...
               + A(:,1)*cott_Q1r2.*S_(2,1,n,xi,eta,Zxi) ...
               + B(:,1)*cott_Q1r2.*T_(2,1,n,eta,Zeta);
    if ESBC && m == M
        indices = logical(r < eps);
        if n == 0
            sigma_pp(:,indices) = repmat(5*(4*a.^2-3*b.^2).*A(:,1)/30,1,sum(indices));
        elseif n == 2
            sigma_pp(:,indices) = repmat((4*a.^2.*A(:,1) - 12*b.^2.*B(:,1))/30,1,sum(indices));
        else
            sigma_pp(:,indices) = 0;
        end
    else
        sigma_pp = sigma_pp + A(:,2)*Q0r2.*S_(6,2,n,xi,eta,Zxi) ...
                            + B(:,2)*Q0r2.*T_(6,2,n,eta,Zeta) ...
                            + A(:,2)*cott_Q1r2.*S_(2,2,n,xi,eta,Zxi) ...
                            + B(:,2)*cott_Q1r2.*T_(2,2,n,eta,Zeta);
    end
	sigma_pp = 2*G*sigma_pp;
    
	solid.sigma_pp = sigma_pp;
end
if options.calc_sigma_rt
    Q1r2 = Q1./r2;
    sigma_rt =   A(:,1)*Q1r2.*S_(7,1,n,xi,eta,Zxi) ...
               + B(:,1)*Q1r2.*T_(7,1,n,eta,Zeta);
    if ESBC && m == M
        sigma_rt(:,logical(r < eps)) = 0;
    else
        sigma_rt = sigma_rt + A(:,2)*Q1r2.*S_(7,2,n,xi,eta,Zxi) ...
                            + B(:,2)*Q1r2.*T_(7,2,n,eta,Zeta);
    end
	sigma_rt = 2*G*sigma_rt;
    
	solid.sigma_rt = sigma_rt;
end
if options.calc_errorsNavier
    Q1sr2 = Q1s./r2;
    sigma_rt =   A(:,1)*Q1sr2.*S_(7,1,n,xi,eta,Zxi) ...
               + B(:,1)*Q1sr2.*T_(7,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        sigma_rt = sigma_rt + A(:,2)*Q1sr2.*S_(7,2,n,xi,eta,Zxi) ...
                            + B(:,2)*Q1sr2.*T_(7,2,n,eta,Zeta);
    end
	sigma_rt = 2*G*sigma_rt;
    
    r3 = r.^3;
    Q0r3 = Q0./r3;
    Q1r3 = Q1./r3;
    Q2r3 = Q2./r3;
    dsigma_rr_dr =   A(:,1)*Q0r3.*S_(8,1,n,xi,eta,Zxi) ...
                   + B(:,1)*Q0r3.*T_(8,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        dsigma_rr_dr =  dsigma_rr_dr + A(:,2)*Q0r3.*S_(8,2,n,xi,eta,Zxi) ...
                                     + B(:,2)*Q0r3.*T_(8,2,n,eta,Zeta);
    end
	dsigma_rr_dr = 2*G*dsigma_rr_dr;
      
    dsigma_rt_dr =   A(:,1)*Q1r3.*S_(9,1,n,xi,eta,Zxi) ...
                   + B(:,1)*Q1r3.*T_(9,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        dsigma_rt_dr =  dsigma_rt_dr + A(:,2)*Q1r3.*S_(9,2,n,xi,eta,Zxi) ...
                                     + B(:,2)*Q1r3.*T_(9,2,n,eta,Zeta);
    end
	dsigma_rt_dr = 2*G*dsigma_rt_dr;

    dsigma_rt_dt =   A(:,1)*Q2r3.*S_(7,1,n,xi,eta,Zxi) ...
                   + B(:,1)*Q2r3.*T_(7,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        dsigma_rt_dt = dsigma_rt_dt + A(:,2)*Q2r3.*S_(7,2,n,xi,eta,Zxi) ...
                                    + B(:,2)*Q2r3.*T_(7,2,n,eta,Zeta);
    end
	dsigma_rt_dt = 2*G*dsigma_rt_dt;

    dsigma_diffr =   A(:,1)*Q1r3.*S_(6,1,n,xi,eta,Zxi) ...
                   + B(:,1)*Q1r3.*T_(6,1,n,eta,Zeta) ...
                   + (-n^2-n+1)*A(:,1)*Q1r3.*S_(2,1,n,xi,eta,Zxi) ...
                   + (-n^2-n+1)*B(:,1)*Q1r3.*T_(2,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        dsigma_diffr = dsigma_diffr + A(:,2)*Q1r3.*S_(6,2,n,xi,eta,Zxi) ...
                                    + B(:,2)*Q1r3.*T_(6,2,n,eta,Zeta) ...
                                    + (-n^2-n+1)*A(:,2)*Q1r3.*S_(2,2,n,xi,eta,Zxi) ...
                                    + (-n^2-n+1)*B(:,2)*Q1r3.*T_(2,2,n,eta,Zeta);
    end
	dsigma_diffr = 2*G*dsigma_diffr;
    
    R = repmat(r,size(xi,1),1);
    Theta = repmat(theta,size(xi,1),1);
    solid.navier1 = dsigma_rr_dr + dsigma_rt_dt + 1./R.*(2*sigma_rr-sigma_tt-sigma_pp+sigma_rt.*cos(Theta));
    solid.navier2 = dsigma_rt_dr + dsigma_diffr + 3./R.*sigma_rt.*sin(Theta);
    if ESBC && m == M
        indices = logical(r < eps);
        if n == 1
            omega = options.omega(options.indices);
            rho_s = options.rho_s(end);
            solid.navier1(:,indices) = repmat(-omega.^2*rho_s.*(a.*A(:,1) - 2*b.*B(:,1))/3,1,sum(indices));
        else
            solid.navier1(:,indices) = 0;
        end
        solid.navier2(:,indices) = 0;
    end
end
