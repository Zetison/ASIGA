function task = getAnalyticSolutions(task)
noDomains = numel(task.varCol);
applyLoad = task.misc.applyLoad;
splitExteriorFields = strcmp(applyLoad,'planeWave') || strcmp(applyLoad,'pointCharge') || strcmp(applyLoad,'radialPulsation');
task.splitExteriorFields = splitExteriorFields;

if ~splitExteriorFields && ~strcmp(task.misc.BC,'NBC')
    error('This exact solution requires BC = ''NBC''')
end    

switch task.misc.applyLoad
    case {'planeWave','radialPulsation'}
        alpha_s = task.ffp.alpha_s;
        beta_s = task.ffp.beta_s;
        d_vec = -[cos(beta_s(1))*cos(alpha_s);
                  cos(beta_s(1))*sin(alpha_s);
                  sin(beta_s(1))*ones(1,length(alpha_s))];
        task.d_vec = d_vec;
    case 'pointCharge'
        alpha_s = task.ffp.alpha_s;
        beta_s = task.ffp.beta_s;
        d_vec =  [cos(beta_s(1))*cos(alpha_s);
                  cos(beta_s(1))*sin(alpha_s);
                  sin(beta_s(1))*ones(1,length(alpha_s))];
        task.d_vec = d_vec;
end
layer = extract_e3Dss_data(task);
if task.analyticSolutionExist
    task.analyticFunctions = @(X) analytic(X,layer);
end
task.p_inc_ = @(X) analytic({X},layer,NaN,'p_inc',1);
task.dp_inc_ = @(X,n) analytic({X},layer,n,'dp_inc',1);
if splitExteriorFields
    task.dp_incdx_ = @(X) analytic({X},layer,NaN,'dp_incdx',1);
    task.dp_incdy_ = @(X) analytic({X},layer,NaN,'dp_incdy',1);
    task.dp_incdz_ = @(X) analytic({X},layer,NaN,'dp_incdz',1);
end
task.dpdn_ = @(X,n) analytic({X},layer,NaN,'dpdn',1);
for i = 1:noDomains
    switch task.varCol{i}.media
        case 'fluid'
            k = task.misc.omega/task.varCol{i}.c_f;
            Phi_k = @(r) exp(1i*k.*r)./(4*pi*r);
            dPhi_kdny = @(xmy,r,ny) -Phi_k(r)./r.^2.*(1i*k.*r - 1).*sum(xmy.*ny,2);
            task.varCol{i}.Phi_k = Phi_k;
            task.varCol{i}.dPhi_kdny = dPhi_kdny;
            
            task.varCol{i}.p_ = @(X) analytic({X},layer,NaN,'p',i);
        case 'solid'
            task.varCol{i}.u_x_ = @(X) analytic({X},layer,NaN,'u_x',i);
            task.varCol{i}.u_y_ = @(X) analytic({X},layer,NaN,'u_y',i);
            task.varCol{i}.u_z_ = @(X) analytic({X},layer,NaN,'u_z',i);
    end
end
task.p_0_ = @(X) analytic({X},layer,NaN,'p_0',1);
task.p_inc_ROM_ = @(X) p_inc_ROM(X,layer,task.p_inc_);
task.dp_inc_ROM_ = @(X,n) dp_inc_ROM(X,n,layer,task.dp_inc_);



function analyticFunctions = analytic(X,varCol,n,func,m_func)
if iscell(X{1})
    X = X{1};
end
M = numel(X);

if isfield(varCol{1},'P_inc')
    P_inc = varCol{1}.P_inc;
else
    P_inc = 1;
end
k = varCol{1}.omega/varCol{1}.c_f;
applyLoad = varCol{1}.applyLoad;
switch applyLoad
    case 'Safjan'
        r = norm2(X{1});
        theta = acos(X{1}(:,3)./r);
        phi = atan2(X{1}(:,2),X{1}(:,1));
        dr = X{1}./r(:,[1,1,1]);
        dtheta = [cos(theta).*cos(phi), cos(theta).*sin(phi), -sin(theta)]./r(:,[1,1,1]);
        varCol{1}.p = hankel_s(1,k*r,1).*cos(theta);
        varCol{1}.p_0 = -cos(theta);
        varCol{1}.dpdx = k*dhankel_s(1,k.*r,1).*cos(theta).*dr(:,1) - hankel_s(1,k.*r,1).*sin(theta).*dtheta(:,1);
        varCol{1}.dpdy = k*dhankel_s(1,k.*r,1).*cos(theta).*dr(:,2) - hankel_s(1,k.*r,1).*sin(theta).*dtheta(:,2);
        varCol{1}.dpdz = k*dhankel_s(1,k.*r,1).*cos(theta).*dr(:,3) - hankel_s(1,k.*r,1).*sin(theta).*dtheta(:,3);
    case 'SimpsonTorus'
        varCol{1}.p = prod(sin(k*X{1}/sqrt(3)),2);
        varCol{1}.p_0 = NaN;
        varCol{1}.dpdx = k/sqrt(3)*cos(k*X{1}(:,1)/sqrt(3)).*sin(k*X{1}(:,2)/sqrt(3)).*sin(k*X{1}(:,3)/sqrt(3));
        varCol{1}.dpdy = k/sqrt(3)*sin(k*X{1}(:,1)/sqrt(3)).*cos(k*X{1}(:,2)/sqrt(3)).*sin(k*X{1}(:,3)/sqrt(3));
        varCol{1}.dpdz = k/sqrt(3)*sin(k*X{1}(:,1)/sqrt(3)).*sin(k*X{1}(:,2)/sqrt(3)).*cos(k*X{1}(:,3)/sqrt(3));
    case 'pointPulsation'
        C_n = @(n) cos(n-1);
        switch varCol{1}.model
            case 'MS'
                R_i = varCol{1}.R_i;
                y = R_i*[1,1,1]/4;
            case 'Barrel'
                R_i = varCol{1}.R;
                y = R_i*[1,1,1]/4;
            case {'BC','BCA'}
                L = varCol{1}.L;
                b = varCol{1}.b;
                x_s = varCol{1}.x_s;
                l_ls = varCol{1}.l_ls;
                h_s = varCol{1}.h_s;
                c = varCol{1}.c;
                arr = b:-b:-(L+2*b);
%                 y = [[x_s-l_ls*0.3,0,c]; [arr; zeros(2,numel(arr))].'];
                y = [arr; zeros(2,numel(arr))].';
%                 y(1) = -100;
            case 'S1_P'
                y = [0,0,0];
%                 y = [1,1,1]/4;
            case 'S1_P2'
                R_i = varCol{1}.R;
                y = R_i/4*[  1,1,1;
                            -1,1,1;
                            1,-1,1;
                            -1,-1,1;
                            1,1,-1;
                            -1,1,-1;
                            1,-1,-1;
                            -1,-1,-1];
            case 'Cube_P'
                a = varCol{1}.a;
                Xarr = linspace(-1,1,3)*a/4;
                Yarr = linspace(-1,1,3)*a/4;
                Zarr = linspace(-1,1,3)*a/4;
                [XX,YY,ZZ] = ndgrid(Xarr,Yarr,Zarr);
                y = [XX(:),YY(:), ZZ(:)];
            otherwise
                y = [0,0,0];
        end
        if ~all(y == 0) && varCol{1}.useROM
            error('Not implemented')
        end
        Phi_k = @(r) exp(1i*k.*r)./(4*pi*r);
        p = zeros(size(X{1},1),numel(k));
        dpdx = zeros(size(X{1},1),numel(k));
        dpdy = zeros(size(X{1},1),numel(k));
        dpdz = zeros(size(X{1},1),numel(k));
        for i = 1:size(y,1)
            xms = X{1}-y(i,:);
            R = norm2(xms);
            p_i = C_n(i)*Phi_k(R);
            p = p + p_i;
            dpdx = dpdx + p_i.*(1i*k - 1./R)./R.*xms(:,1);
            dpdy = dpdy + p_i.*(1i*k - 1./R)./R.*xms(:,2);
            dpdz = dpdz + p_i.*(1i*k - 1./R)./R.*xms(:,3);
        end
        p_0 = zeros(size(X{1},1),1);
        for i = 1:size(y,1)
            p_0 = p_0 + C_n(i)/(4*pi)*exp(-1i*k.*dot3(X{1}./repmat(norm2(X{1}),1,3),y(i,:).'));
        end
        varCol{1}.p = p;
        varCol{1}.p_0 = p_0;
        varCol{1}.dpdx = dpdx;
        varCol{1}.dpdy = dpdy;
        varCol{1}.dpdz = dpdz;
    case {'planeWave','radialPulsation','pointCharge'}
        d_vec = varCol{1}.d_vec;
        if strcmp(applyLoad,'pointCharge')
            r_s = varCol{1}.r_s;
            options.r_s = r_s;
        end
        
        isSphericalShell = varCol{1}.isSphericalShell;
        N_max = varCol{1}.N_max;
        if ~(nargin < 3) && (strcmp(func,'p_inc') || strcmp(func,'dp_inc') || strcmp(func,'dp_incdx') || strcmp(func,'dp_incdy') || strcmp(func,'dp_incdz')) && isinf(N_max)
            if strcmp(applyLoad,'planeWave')
                varCol{1}.p_inc = P_inc*exp(1i*(X{1}*d_vec).*k);
                if nargin > 2 && ~any(isnan(n(:)))
                    varCol{1}.dp_inc = 1i*(n*d_vec).*k.*varCol{1}.p_inc;
                end
                varCol{1}.dp_incdx = 1i*k.*varCol{1}.p_inc.*d_vec(1,:);
                varCol{1}.dp_incdy = 1i*k.*varCol{1}.p_inc.*d_vec(2,:);
                varCol{1}.dp_incdz = 1i*k.*varCol{1}.p_inc.*d_vec(3,:);
            elseif strcmp(applyLoad,'radialPulsation')
                R_i = varCol{1}.R_i;
                varCol{1}.p_inc = P_inc*R_i.*exp(-1i*k.*(norm2(X{1})-R_i))./norm2(X{1});
                R = norm2(X{1});
                if nargin > 2 && ~any(isnan(n(:)))
                    varCol{1}.dp_inc = -varCol{1}.p_inc.*(1i*k+1./R).*sum(X{1}.*n,2)./R;
                end
                varCol{1}.dp_incdx = -varCol{1}.p_inc.*(1i*k+1./R).*X{1}(:,1)./R;
                varCol{1}.dp_incdy = -varCol{1}.p_inc.*(1i*k+1./R).*X{1}(:,2)./R;
                varCol{1}.dp_incdz = -varCol{1}.p_inc.*(1i*k+1./R).*X{1}(:,3)./R;
            elseif strcmp(applyLoad,'pointCharge')
                x_s = r_s*d_vec.';
                xms = X{1}-x_s;
                R = norm2(xms);
                varCol{1}.p_inc = P_inc*r_s.*exp(1i*k.*R)./R;
                if nargin > 2 && ~any(isnan(n(:)))
                    varCol{1}.dp_inc = varCol{1}.p_inc.*(1i*k-1./R).*sum(xms.*n,2)./R;
                end
                varCol{1}.dp_incdx = varCol{1}.p_inc.*(1i*k-1./R).*xms(:,1)./R;
                varCol{1}.dp_incdy = varCol{1}.p_inc.*(1i*k-1./R).*xms(:,2)./R;
                varCol{1}.dp_incdz = varCol{1}.p_inc.*(1i*k-1./R).*xms(:,3)./R;
            end
        else
            options.BC = varCol{1}.BC;
            options.d_vec = d_vec;
            options.N_max = N_max;
            options.P_inc = P_inc;
            options.omega = varCol{1}.omega;
            options.applyLoad = applyLoad;
            options.Display = 'none';
            if nargin < 3
                calc_p_0 = false;
            else
                calc_p_0 = strcmp(func,'p_0');
            end
            varCol{1}.calc_p_0 = calc_p_0 && isSphericalShell && ~isempty(X{1});      % Toggle calculation of the far field pattern
            for m = 1:M
                varCol{m}.X = X{m};

                % Parameters in layer m for options{i}.media = 'fluid'
                switch varCol{m}.media
                    case 'fluid'
                        varCol{m}.calc_p       	 = isSphericalShell.*~calc_p_0;      % Toggle calculation of the scattered pressure
                        varCol{m}.calc_dp      	 = isSphericalShell.*true(1,3).*~calc_p_0; % Toggle calculation of the three components of the gradient of the pressure
                        varCol{m}.calc_p_laplace = isSphericalShell.*~calc_p_0;      % Toggle calculation of the Laplace operator of the scattered pressure fields
                        varCol{m}.calc_errHelm	 = isSphericalShell.*~calc_p_0;      % Toggle calculation of the errors for the Helmholtz equation
                        varCol{m}.calc_p_inc     = true.*~calc_p_0;      % Toggle calculation of the incident pressure
                        varCol{m}.calc_dp_inc    = true(1,3).*~calc_p_0; % Toggle calculation of the three components of the gradient of the incident pressure
                    case 'solid'
                        % Parameters in layer m for options{i}.media = 'solid' or 'viscoelastic'
                        varCol{m}.calc_u       = isSphericalShell.*true(1,3).*~calc_p_0; % Toggle calculation of the three components of the displacement
                        varCol{m}.calc_du      = isSphericalShell.*true(3,3).*~calc_p_0; % Toggle calculation of the three cartesian derivatives of the three components of the displacement [du_xdx du_xdy du_xdz; 
                                                                              %                                                                                                    du_ydx du_ydy du_ydz; 
                                                                              %                                                                                                    du_zdx du_zdy du_zdz]
                        varCol{m}.calc_sigma   = isSphericalShell.*true(1,6).*~calc_p_0; % Toggle calculation of the six components of the stress field (cartesian coordinates) [sigma_xx sigma_yy sigma_zz sigma_yz sigma_xz sigma_xy]
                end
            end
            varCol = e3Dss(varCol, options);
        end
    case 'Cartesian'
        k = varCol{1}.omega/varCol{1}.c_f;


        Acoeff = [1, 1i+1];
        Bcoeff = [1-1i, 1i];
        Ccoeff = [0.5, 1i+1];
        Dcoeff = [1-1i, 1i];

        Ecoeff = [1-1i, 1i;
                  1,    1];

        Fcoeff = [1+1i, 1;
                  1i,    -1];
        [p,dpdx,dpdy,dpdz] = general3DSolutionHelmholtz(X{1}, k, Acoeff,Bcoeff,Ccoeff,Dcoeff,Ecoeff,Fcoeff);
        varCol{1}.p = p;
        varCol{1}.dpdx = dpdx;
        varCol{1}.dpdy = dpdy;
        varCol{1}.dpdz = dpdz;
end

if nargin > 2 && ~any(isnan(n(:)))
    if isfield(varCol{1},'dpdx')
        varCol{1}.dpdn = varCol{1}.dpdx.*n(:,1)+varCol{1}.dpdy.*n(:,2)+varCol{1}.dpdz.*n(:,3);
    end
    if ~isfield(varCol{1},'dp_inc')
        if varCol{1}.splitExteriorFields
            varCol{1}.dp_inc = varCol{1}.dp_incdx.*n(:,1)+varCol{1}.dp_incdy.*n(:,2)+varCol{1}.dp_incdz.*n(:,3);
        else
            varCol{1}.dp_inc = -varCol{1}.dpdn;
        end
    end
end
if nargin < 3
    analyticFunctions = varCol; 
else
    analyticFunctions = varCol{m_func}.(func);
end


function [p,dpdx,dpdy,dpdz] = general3DSolutionHelmholtz(X, k, A, B, C, D, E, F)
% Solution with notation from 
% http://mathworld.wolfram.com/HelmholtzDifferentialEquationCartesianCoordinates.html
p = 0;
x = X(:,1);
y = Y(:,1);
z = Z(:,1);
dpdx = zeros(numel(x),numel(k));
dpdy = zeros(numel(x),numel(k));
dpdz = zeros(numel(x),numel(k));
[L, M] = size(E);

for l = 1:L
    for m = 1:M
        lambda = sqrt(k.^2+l^2+m^2);
        p = p + (A(l)*exp(l*x) + B(l)*exp(-l*x)).*(C(m)*exp(m*y) + D(m)*exp(-m*y)).*(E(l,m)*exp(-1i*lambda.*z) + F(l,m)*exp(1i*lambda.*z));
        dpdx = dpdx + (A(l)*l*exp(l*x) - B(l)*l*exp(-l*x)).*(C(m)*exp(m*y) + D(m)*exp(-m*y)).*(E(l,m)*exp(-1i*lambda.*z) + F(l,m)*exp(1i*lambda.*z));
        dpdy = dpdy + (A(l)*exp(l*x) + B(l)*exp(-l*x)).*(C(m)*m*exp(m*y) - D(m)*m*exp(-m*y)).*(E(l,m)*exp(-1i*lambda.*z) + F(l,m)*exp(1i*lambda.*z));
        dpdz = dpdz + (A(l)*exp(l*x) + B(l)*exp(-l*x)).*(C(m)*exp(m*y) + D(m)*exp(-m*y)).*(E(l,m)*-1i*lambda.*exp(-1i*lambda.*z) + F(l,m)*1i*lambda.*exp(1i*lambda.*z));
    end
end



function dp_inc_ROM = dp_inc_ROM(X,n,varCol,dp_inc_)

m = 0:(varCol{1}.noRHSs-1);
switch varCol{1}.applyLoad
    case 'planeWave'
        d_vecX = X*varCol{1}.d_vec;
    case 'pointPulsation'
        d_vecX = norm2(X);
    case 'pointCharge'
        x_s = varCol{1}.r_s*varCol{1}.d_vec.';
        d_vecX = norm2(X-x_s);
    otherwise
        error('Not implemented')
end
temp = zeros(numel(d_vecX),varCol{1}.noRHSs);
c_f = varCol{1}.c_f;
k = varCol{1}.omega/c_f;
temp(:,2:end) = (1i./d_vecX)*m(2:end)/k;
dp_inc_ROM = dp_inc_(X,n).*(1i*d_vecX/c_f).^m.*(1-temp);

function p_inc_ROM = p_inc_ROM(X,varCol,p_inc_)

m = 0:(varCol{1}.noRHSs-1);
switch varCol{1}.applyLoad
    case 'planeWave'
        d_vecX = X*varCol{1}.d_vec;
    case 'pointPulsation'
        d_vecX = norm2(X);
    case 'pointCharge'
        x_s = varCol{1}.r_s*varCol{1}.d_vec.';
        d_vecX = norm2(X-x_s);
    otherwise
        error('Not implemented')
end
c_f = varCol{1}.c_f;
p_inc_ROM = p_inc_(X).*(1i*d_vecX/c_f).^m;





