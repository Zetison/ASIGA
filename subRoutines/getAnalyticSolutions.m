function varCol = getAnalyticSolutions(varCol)
k = varCol{1}.k;
omega = varCol{1}.omega;
if ~isfield(varCol{1},'P_inc')
    P_inc = 1;
else
    P_inc = varCol{1}.P_inc;
end
BC = varCol{1}.BC;
applyLoad = varCol{1}.applyLoad;
switch applyLoad
    case {'pointPulsation','SimpsonTorus','Safjan'}        
        p_inc = @(v) zeros(size(v,1),1);
        gp_inc = @(v) zeros(size(v,1),1);
        dpdn = @(X,n) analytic({X},varCol,n,'dpdn',1);
        dp_inc = @(v,n) -dpdn(v,n);
        if ~strcmp(BC,'NBC')
            error('This exact solution is not a solution to coupled problems')
        end        
    case {'planeWave','radialPulsation'}
        alpha_s = varCol{1}.alpha_s;
        beta_s = varCol{1}.beta_s;
        d_vec = -[cos(beta_s(1))*cos(alpha_s);
                  cos(beta_s(1))*sin(alpha_s);
                  sin(beta_s(1))*ones(1,length(alpha_s))];
        varCol{1}.d_vec = d_vec;
        N_max = varCol{1}.N_max;
        e3Dss_options.BC = BC;
        e3Dss_options.d_vec = d_vec;
        e3Dss_options.N_max = N_max;
        e3Dss_options.P_inc = P_inc;
        e3Dss_options.omega = omega;
        e3Dss_options.applyLoad = applyLoad;
        e3Dss_options.Display = 'none';
        varCol{1}.e3Dss_options = e3Dss_options;
        if strcmp(applyLoad,'planeWave')
            if isinf(N_max)
                p_inc = @(v) P_inc*exp(1i*(v*d_vec)*k);
                gp_inc = @(v) 1i*elementProd(p_inc(v),d_vec.')*k;
                dp_inc = @(v,n) 1i*(n*d_vec)*k.*p_inc(v);
            else
                p_inc = @(v) P_inc*exp(1i*(v*d_vec)*k);
                gp_inc = @(X) -analytic({X},varCol,NaN,'dp',1);
                dp_inc = @(X,n) -analytic({X},varCol,n,'dpdn',1);
            end
        elseif strcmp(applyLoad,'radialPulsation')
            p_inc = @(v) P_inc*R_o(1).*exp(-1i*k*(norm2(v)-R_o(1)))./norm2(v);
            gp_inc = @(v)   -p_inc(v).*(1./norm2(v)+1i*k).*v./norm2(v);
            dp_inc = @(v,n) -p_inc(v).*(1./norm2(v)+1i*k).*sum(v.*n,2)./norm2(v);
        elseif strcmp(applyLoad,'pointCharge')
            r_s = varCol{1}.r_s;
            x_s = d_vec*r_s;
            p_inc = @(v) P_inc*r_s.*exp(1i*k*norm2(v-x_s))./norm2(v-x_s);
            gp_inc = @(v)   p_inc(v).*(-1./norm2(v-x_s)+1i*k).*v./norm2(v-x_s);
            dp_inc = @(v,n) p_inc(v).*(-1./norm2(v-x_s)+1i*k).*sum((v-x_s).*n,2)./norm2(v-x_s);
        end
        dpdn = @(v,n) zeros(size(v,1),1);
        varCol{1}.d_vec = d_vec;
end
Phi_k = @(r) exp(1i*k*r)./(4*pi*r);
dPhi_kdny = @(xmy,r,ny) -Phi_k(r)./r.^2.*(1i*k*r - 1).*sum(xmy.*ny,2);
varCol{1}.Phi_k = Phi_k;
varCol{1}.dPhi_kdny = dPhi_kdny;
if varCol{1}.analyticSolutionExist
    varCol{1}.analyticFunctions = @(X) analytic(X,varCol);
end
noRHSs = varCol{1}.noRHSs;
c_f = varCol{1}.c_f;
varCol{1}.p_inc = p_inc;
varCol{1}.dp_inc = dp_inc;
varCol{1}.p_inc_ROM = @(X) p_inc_ROM(X,c_f,d_vec,noRHSs,p_inc);
varCol{1}.dp_inc_ROM = @(X,n) dp_inc_ROM(X,n,omega,c_f,d_vec,noRHSs,p_inc);
varCol{1}.gp_inc = gp_inc;
varCol{1}.dpdn = dpdn;
varCol{1}.p = @(X) analytic({X},varCol,NaN,'p',1);
varCol{1}.p_0 = @(X) analytic({X},varCol,NaN,'p_0',1);



function analyticFunctions = analytic(X,layers,n,func,m_func)
M = numel(X);
switch layers{1}.applyLoad
    case 'Safjan'
        k = layers{1}.k;
        r = norm2(X{1});
        theta = acos(X{1}(:,3)./r);
        phi = atan2(X{1}(:,2),X{1}(:,1));
        dr = X{1}./r(:,[1,1,1]);
        dtheta = [cos(theta).*cos(phi), cos(theta).*sin(phi), -sin(theta)]./r(:,[1,1,1]);
        layers{1}.p = hankel_s(1,k*r,1).*cos(theta);
        layers{1}.p_0 = -cos(theta);
        layers{1}.dp = k*repmat(dhankel_s(1,k*r,1).*cos(theta),1,3).*dr - repmat(hankel_s(1,k*r,1).*sin(theta),1,3).*dtheta;
    case 'SimpsonTorus'
        k = layers{1}.k;
        layers{1}.p = prod(sin(k*X{1}/sqrt(3)),2);
        layers{1}.p_0 = NaN;
        layers{1}.dp = k/sqrt(3)*[cos(k*X{1}(:,1)/sqrt(3)).*sin(k*X{1}(:,2)/sqrt(3)).*sin(k*X{1}(:,3)/sqrt(3)), ...
                                    sin(k*X{1}(:,1)/sqrt(3)).*cos(k*X{1}(:,2)/sqrt(3)).*sin(k*X{1}(:,3)/sqrt(3)), ...
                                    sin(k*X{1}(:,1)/sqrt(3)).*sin(k*X{1}(:,2)/sqrt(3)).*cos(k*X{1}(:,3)/sqrt(3))];
    case 'pointPulsation'
        k = layers{1}.k;
        C_n = @(n) cos(n-1);
        switch layers{1}.model
            case 'MS'
                R = layers{1}.R;
                y = R*[1,1,1]/4;
            case 'Barrel'
                R = layers{1}.R;
                y = R*[1,1,1]/4;
            case {'BC','BCA'}
                L = layers{1}.L;
                b = layers{1}.b;
                x_s = layers{1}.x_s;
                l_ls = layers{1}.l_ls;
                h_s = layers{1}.h_s;
                c = layers{1}.c;
                arr = b:-b:-(L+2*b);
%                 y = [[x_s-l_ls*0.3,0,c]; [arr; zeros(2,numel(arr))].'];
                y = [arr; zeros(2,numel(arr))].';
%                 y(1) = -100;
            case 'S1_P'
                y = [0,0,0];
%                 y = [1,1,1]/4;
            case 'S1_P2'
                R = layers{1}.R;
                y = R/4*[1,1,1;
                        -1,1,1;
                        1,-1,1;
                        -1,-1,1;
                        1,1,-1;
                        -1,1,-1;
                        1,-1,-1;
                        -1,-1,-1];
            case 'Cube_P'
                a = layers{1}.a;
                Xarr = linspace(-1,1,3)*a/4;
                Yarr = linspace(-1,1,3)*a/4;
                Zarr = linspace(-1,1,3)*a/4;
                [XX,YY,ZZ] = ndgrid(Xarr,Yarr,Zarr);
                y = [XX(:),YY(:), ZZ(:)];
            otherwise
                y = [0,0,0];
        end
        Phi_k = @(r) exp(1i*k(1)*r)./(4*pi*r);
        p = zeros(size(X{1},1),1);
        dp = zeros(size(X{1},1),3);
        for i = 1:size(y,1)
            xms = @(v) [v(:,1)-y(i,1),v(:,2)-y(i,2),v(:,3)-y(i,3)];
            R = @(v) norm2(xms(v));
            p_i = C_n(i)*Phi_k(R(X{1}));
            p = p + p_i;
            dp = dp + elementProd(p_i.*(1i*k - 1./R(X{1}))./R(X{1}), xms(X{1}));
        end
        p_0 = zeros(size(X{1},1),1);
        for i = 1:size(y,1)
            p_0 = p_0 + C_n(i)/(4*pi)*exp(-1i*k*dot3(X{1}./repmat(norm2(X{1}),1,3),y(i,:).'));
        end
        layers{1}.p = p;
        layers{1}.p_0 = p_0;
        layers{1}.dp = dp;
        layers{1}.dpdx = dp(:,1);
        layers{1}.dpdy = dp(:,2);
        layers{1}.dpdz = dp(:,3);
        
    case {'planeWave','radialPulsation'}
        layers{1}.calc_p_0       = true;      % Toggle calculation of the far field pattern
        for m = 1:M
            layers{m}.X = X{m};

            % Parameters in layer m for options{i}.media = 'fluid'
            layers{m}.calc_p       	= true;      % Toggle calculation of the scattered pressure
            layers{m}.calc_dp      	= true(1,3); % Toggle calculation of the three components of the gradient of the pressure
            layers{m}.calc_p_laplace	= true;      % Toggle calculation of the Laplace operator of the scattered pressure fields
            layers{m}.calc_errHelm	= true;      % Toggle calculation of the errors for the Helmholtz equation

            % Parameters in layer m for options{i}.media = 'solid' or 'viscoelastic'
            layers{m}.calc_u       = true(1,3); % Toggle calculation of the three components of the displacement
            layers{m}.calc_du      = true(3,3); % Toggle calculation of the three cartesian derivatives of the three components of the displacement [du_xdx du_xdy du_xdz; 
                                                %                                                                                                    du_ydx du_ydy du_ydz; 
                                                %                                                                                                    du_zdx du_zdy du_zdz]
            layers{m}.calc_sigma   = true(1,6); % Toggle calculation of the six components of the stress field (cartesian coordinates) [sigma_xx sigma_yy sigma_zz sigma_yz sigma_xz sigma_xy]
        end
        layers = e3Dss(layers, layers{1}.e3Dss_options);
        if nargin > 2
            layers{1}.dp = [layers{m}.dpdx,layers{m}.dpdy,layers{m}.dpdz];
        end
    case 'Cartesian'
        k = layers{1}.c_f/layers{m}.e3Dss_options.omega;


        Acoeff = [1, 1i+1];
        Bcoeff = [1-1i, 1i];
        Ccoeff = [0.5, 1i+1];
        Dcoeff = [1-1i, 1i];

        Ecoeff = [1-1i, 1i;
                  1,    1];

        Fcoeff = [1+1i, 1;
                  1i,    -1];
        [p,dp] = general3DSolutionHelmholtz(X{m}, k, Acoeff,Bcoeff,Ccoeff,Dcoeff,Ecoeff,Fcoeff);
        layers{1}.p = p;
        layers{1}.dp = dp;
        layers{1}.dpdx = dp(:,1);
        layers{1}.dpdy = dp(:,2);
        layers{1}.dpdz = dp(:,3);
end
if nargin > 2 && ~any(isnan(n(:)))
    for m = 1:M
        layers{m}.dpdn = sum(layers{m}.dp.*n,2);
    end
end
if nargin < 3
    analyticFunctions = layers; 
else
    analyticFunctions = layers{m_func}.(func);
end


function [p,dp] = general3DSolutionHelmholtz(X, k, A, B, C, D, E, F)
% Solution with notation from 
% http://mathworld.wolfram.com/HelmholtzDifferentialEquationCartesianCoordinates.html
p = 0;
dp = zeros(3,1);
x = X(:,1);
y = Y(:,1);
z = Z(:,1);
[L, M] = size(E);

for l = 1:L
    for m = 1:M
        lambda = sqrt(k^2+l^2+m^2);
        p = p + (A(l)*exp(l*x) + B(l)*exp(-l*x)).*(C(m)*exp(m*y) + D(m)*exp(-m*y)).*(E(l,m)*exp(-1i*lambda*z) + F(l,m)*exp(1i*lambda*z));
        dp(1) = dp(1) + (A(l)*l*exp(l*x) - B(l)*l*exp(-l*x)).*(C(m)*exp(m*y) + D(m)*exp(-m*y)).*(E(l,m)*exp(-1i*lambda*z) + F(l,m)*exp(1i*lambda*z));
        dp(2) = dp(2) + (A(l)*exp(l*x) + B(l)*exp(-l*x)).*(C(m)*m*exp(m*y) - D(m)*m*exp(-m*y)).*(E(l,m)*exp(-1i*lambda*z) + F(l,m)*exp(1i*lambda*z));
        dp(3) = dp(3) + (A(l)*exp(l*x) + B(l)*exp(-l*x)).*(C(m)*exp(m*y) + D(m)*exp(-m*y)).*(E(l,m)*-1i*lambda*exp(-1i*lambda*z) + F(l,m)*1i*lambda*exp(1i*lambda*z));
    end
end



function dp_inc_ROM = dp_inc_ROM(X,n,omega,c_f,d_vec,noRHSs,p_inc)

m = 0:(noRHSs-1);
d_vecX = X*d_vec;
temp = zeros(numel(d_vecX),noRHSs);
temp(:,2:end) = (1./d_vecX)*m(2:end);
dp_inc_ROM = (n*d_vec).*p_inc(X).*(1i*d_vecX/c_f).^m.*(temp+1i*omega/c_f);

function p_inc_ROM = p_inc_ROM(X,c_f,d_vec,noRHSs,p_inc)

m = 0:(noRHSs-1);
d_vecX = X*d_vec;
p_inc_ROM = p_inc(X).*(1i*d_vecX/c_f).^m;





