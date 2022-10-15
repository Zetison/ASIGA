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
p_inc_ = @(X) analytic({X},layer,NaN,'p_inc',1);
dp_incdn_ = @(X,n) analytic({X},layer,n,'dp_incdn',1);
task.p_inc_ = p_inc_;
task.dp_incdn_ = dp_incdn_;
if splitExteriorFields
    task.dp_incdx_ = @(X) analytic({X},layer,NaN,'dp_incdx',1);
    task.dp_incdy_ = @(X) analytic({X},layer,NaN,'dp_incdy',1);
    task.dp_incdz_ = @(X) analytic({X},layer,NaN,'dp_incdz',1);
end
task.dpdn_ = @(X,n) analytic({X},layer,n,'dpdn',1);
for i = 1:noDomains
    switch task.varCol{i}.media
        case 'fluid'
            k = task.misc.omega/task.varCol{i}.c_f;
            Phi_k = @(r) exp(1i*k.*r)./(4*pi*r);
            dPhi_kdny = @(xmy,r,ny) -Phi_k(r)./r.^2.*(1i*k.*r - 1).*sum(xmy.*ny,2);
            dPhi_kdnx = @(xmy,r,ny)  Phi_k(r)./r.^2.*(1i*k.*r - 1).*sum(xmy.*ny,2);
            task.varCol{i}.Phi_k = Phi_k;
            task.varCol{i}.dPhi_kdny = dPhi_kdny;
            task.varCol{i}.dPhi_kdnx = dPhi_kdnx;
            
            task.varCol{i}.p_ = @(X) analytic({X},layer,NaN,'p',i);
        case 'solid'
            task.varCol{i}.u_ = @(X) analytic({X},layer,NaN,'u',i);
    end
end
task.p_0_ = @(X) analytic({X},layer,NaN,'p_0',1);
omega = task.misc.omega;
symmetric = task.misc.symmetric;
task.p_inc_ROM_ = @(X) p_inc_ROM(X,layer,p_inc_,omega,symmetric);
task.dp_incdn_ROM_ = @(X,n) dp_incdn_ROM(X,n,layer,dp_incdn_);



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
        varCol{1}.p = sphBessel(1,k.*r,3).*cos(theta);
        varCol{1}.p_0 = -cos(theta)./k;
        varCol{1}.dp = cell(1,3);
        for i = 1:3
            varCol{1}.dp{i} = k.*dSphBessel(1,k.*r,3).*cos(theta).*dr(:,i) - sphBessel(1,k.*r,3).*sin(theta).*dtheta(:,i);
        end
    case 'Safjan10'
        r = norm2(X{1});
        theta = acos(X{1}(:,3)./r);
        phi = atan2(X{1}(:,2),X{1}(:,1));
        dr = X{1}./r(:,[1,1,1]);
        dtheta = [cos(theta).*cos(phi), cos(theta).*sin(phi), -sin(theta)]./r(:,[1,1,1]);
        P = legendre(10,cos(theta)); 
        P = P(1,:).';
        Pnm1 = legendre(9,cos(theta));
        dP = 10*(P.*cos(theta) - Pnm1(1,:).')./sin(theta);
        varCol{1}.p = sphBessel(10,k.*r,3).*P;
        varCol{1}.p_0 = -1i*P./k;
        varCol{1}.dp = cell(1,3);
        for i = 1:3
            varCol{1}.dp{i} = k.*dSphBessel(10,k.*r,3).*P.*dr(:,i) - sphBessel(10,k.*r,3).*dP.*dtheta(:,i);
        end
    case 'SimpsonTorus'
        varCol{1}.p = prod(sin(k*X{1}/sqrt(3)),2);
        varCol{1}.p_0 = NaN;
        
        varCol{1}.dp{1} = k/sqrt(3)*cos(k*X{1}(:,1)/sqrt(3)).*sin(k*X{1}(:,2)/sqrt(3)).*sin(k*X{1}(:,3)/sqrt(3));
        varCol{1}.dp{2} = k/sqrt(3)*sin(k*X{1}(:,1)/sqrt(3)).*cos(k*X{1}(:,2)/sqrt(3)).*sin(k*X{1}(:,3)/sqrt(3));
        varCol{1}.dp{3} = k/sqrt(3)*sin(k*X{1}(:,1)/sqrt(3)).*sin(k*X{1}(:,2)/sqrt(3)).*cos(k*X{1}(:,3)/sqrt(3));
        
    case 'pointPulsation'
        C_n = @(n) cos(n-1);
        switch varCol{1}.model
            case 'MS'
                R = varCol{1}.R;
                y = R*[1,1,1]/4;
            case {'Barrel','Barrel_Sweep'}
                R = varCol{1}.R;
                Xarr = linspace(-1,1,3)*R/4;
                Yarr = linspace(-1,1,3)*R/4;
                Zarr = linspace(-1,1,3)*R/4;
                [XX,YY,ZZ] = ndgrid(Xarr,Yarr,Zarr);
                y = [XX(:),YY(:), ZZ(:)];
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
                R = varCol{1}.R;
                y = R/4*[  1,1,1;
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
%         if ~all(y(:) == 0) && varCol{1}.rom.useROM
%             error('Not implemented')
%         end
        Phi_k = @(r) exp(1i*k.*r)./(4*pi*r);
        p = zeros(size(X{1},1),numel(k));
        dp = cell(1,3);
        for i = 1:3
            dp{i} = zeros(size(X{1},1),numel(k));
        end
        for i = 1:size(y,1)
            xms = X{1}-y(i,:);
            R = norm2(xms);
            p_i = C_n(i)*Phi_k(R);
            p = p + p_i;
            for j = 1:3
                dp{j} = dp{j} + p_i.*(1i*k - 1./R)./R.*xms(:,j);
            end
        end
        p_0 = zeros(size(X{1},1),1);
        for i = 1:size(y,1)
            p_0 = p_0 + C_n(i)/(4*pi)*exp(-1i*k.*dot3(X{1}./repmat(norm2(X{1}),1,3),y(i,:).'));
        end
        varCol{1}.p = p;
        varCol{1}.p_0 = p_0;
        varCol{1}.dp = dp;
    case {'planeWave','radialPulsation','pointCharge'}
        d_vec = varCol{1}.d_vec;
        if strcmp(applyLoad,'pointCharge')
            r_s = varCol{1}.r_s;
            options.r_s = r_s;
        end
        
        isSphericalShell = varCol{1}.isSphericalShell;
        N_max = varCol{1}.N_max;
        if ~(nargin < 3) && (strcmp(func,'p_inc') || strcmp(func,'dp_inc') || strcmp(func,'dp_incdn')) && isinf(N_max)
            if strcmp(applyLoad,'planeWave')
                varCol{1}.p_inc = P_inc*exp(1i*(X{1}*d_vec).*k);
                for i = 1:3
                    varCol{1}.dp_inc{i} = 1i*k.*varCol{1}.p_inc.*d_vec(i,:);
                end
            elseif strcmp(applyLoad,'radialPulsation')
                R = varCol{1}.R;
                varCol{1}.p_inc = P_inc*R.*exp(-1i*k.*(norm2(X{1})-R))./norm2(X{1});
                R = norm2(X{1});
                for i = 1:3
                    varCol{1}.dp_inc{i} = -varCol{1}.p_inc.*(1i*k+1./R).*X{1}(:,i)./R;
                end
            elseif strcmp(applyLoad,'pointCharge')
                x_s = r_s*d_vec.';
                xms = X{1}-x_s;
                R = norm2(xms);
                varCol{1}.p_inc = P_inc*r_s.*exp(1i*k.*R)./R;
                for i = 1:3
                    varCol{1}.dp_inc{i} = varCol{1}.p_inc.*(1i*k-1./R).*xms(:,i)./R;
                end
            end
        else
            if isinf(N_max) && nargin >= 3 && (strcmp(func,'p_inc') || strcmp(func,'dp_incdn'))
                for m = 1:M
                    k = varCol{1}.omega/varCol{1}.c_f;
                    if strcmp(applyLoad,'planeWave')
                        k_vec = k.'*d_vec.';
                        varCol{1}.p_inc = P_inc*exp(1i*k_vec*X{m}.').';
                        for ii = 1:3
                            varCol{m}.dp_inc{ii} = 1i*k_vec(:,ii).'.*varCol{1}.p_inc;
                        end
                    else
                        x_s = r_s*d_vec.';
                        Xxms = X{m}-x_s;
                        nXxms = norm2(Xxms);
                        p_inc = P_inc.*exp(1i*k.*nXxms)./nXxms;
                        varCol{m}.p_inc = p_inc;
                        for ii = 1:3
                            varCol{m}.dp_inc{ii} = p_inc.*(1i*k - 1./nXxms).*Xxms(:,ii)./nXxms;
                        end
                    end
                end
            else
                options.BC = varCol{1}.BC;
                options.d_vec = d_vec;
                options.N_max = N_max;
                options.P_inc = P_inc;
                options.omega = varCol{1}.omega;
                options.applyLoad = applyLoad;
                options.Display = 'none';
                options.p_inc_fromSeries = ~isinf(N_max);
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
        varCol{1}.dp{1} = dpdx;
        varCol{1}.dp{2} = dpdy;
        varCol{1}.dp{3} = dpdz;
end

if nargin > 2 && ~any(isnan(n(:)))
    if isfield(varCol{1},'dp')
        varCol{1}.dpdn = varCol{1}.dp{1}.*n(:,1)+varCol{1}.dp{2}.*n(:,2)+varCol{1}.dp{3}.*n(:,3);
    end
    if ~isfield(varCol{1},'dp_incdn')
        if varCol{1}.splitExteriorFields
            varCol{1}.dp_incdn = varCol{1}.dp_inc{1}.*n(:,1)+varCol{1}.dp_inc{2}.*n(:,2)+varCol{1}.dp_inc{3}.*n(:,3);
%             varCol{1}.dp_inc = -varCol{1}.dpdn;
        else
            varCol{1}.dp_incdn = -varCol{1}.dpdn;
        end
    end
end
if nargin < 3
    analyticFunctions = varCol; 
else
    switch func
        case 'dp_incdx'
            varCol{m_func}.(func) = varCol{1}.dp_inc{1};
        case 'dp_incdy'
            varCol{m_func}.(func) = varCol{1}.dp_inc{2};
        case 'dp_incdz'
            varCol{m_func}.(func) = varCol{1}.dp_inc{3};
    end
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



function dp_inc_ROM = dp_incdn_ROM(X,n,varCol,dp_incdn_)

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
dp_inc_ROM = dp_incdn_(X,n).*(1i*d_vecX/c_f).^m.*(1-temp);

function p_inc_ROM = p_inc_ROM(X,varCol,p_inc_,omega,symmetric)

if symmetric
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
    mm1 = m - 1;
    mm1(1) = 0;
    mm2 = m - 2;
    mm2(1) = 0;
    if numel(mm2) > 1
        mm2(2) = 0;
    end
    p_inc_ROM = p_inc_(X).*(omega.^2.*(1i*d_vecX/c_f).^m + 2*m.*omega.*(1i*d_vecX/c_f).^mm1 + m.*(m-1).*(1i*d_vecX/c_f).^mm2);
else
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
end
    
    

