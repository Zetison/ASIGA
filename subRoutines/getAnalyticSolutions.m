
Phi_k = @(r) exp(1i*k_1*r)./(4*pi*r);
dPhi_kdny = @(xmy,r,ny) -Phi_k(r)./r.^2.*(1i*k_1*r - 1).*sum(xmy.*ny,2);
dPhi_kdnx = @(xmy,r,ny) Phi_k(r)./r.^2.*(1i*k_1*r - 1).*sum(xmy.*ny,2);
varCol{1}.Phi_k = Phi_k;
varCol{1}.dPhi_kdny = dPhi_kdny;

switch applyLoad
    case 'SimpsonTorus'
        analytic = @(v) prod(sin(k_1*v/sqrt(3)),2);
        gAnalytic = @(v) k_1/sqrt(3)*[cos(k_1*v(:,1)/sqrt(3)).*sin(k_1*v(:,2)/sqrt(3)).*sin(k_1*v(:,3)/sqrt(3)), ...
                                    sin(k_1*v(:,1)/sqrt(3)).*cos(k_1*v(:,2)/sqrt(3)).*sin(k_1*v(:,3)/sqrt(3)), ...
                                    sin(k_1*v(:,1)/sqrt(3)).*sin(k_1*v(:,2)/sqrt(3)).*cos(k_1*v(:,3)/sqrt(3))];
        dpdn = @(v,n) sum(gAnalytic(v).*n,2);
        
        varCol{1}.farField = @(v) NaN;
        varCol{1}.analytic = analytic;
        varCol{1}.gAnalytic = gAnalytic;
        
        p_inc = @(v) zeros(size(v,1),1);
        gp_inc = @(v) zeros(size(v,1),1);
        dp_inc = @(v,n) -dpdn(v,n);
        
        e3Dss_options = struct('omega', omega, ...
                         'P_inc', P_inc, ...
                         'rho_f', rho_f, ...
                         'c_f', c_f);
    case 'radialPulsation'
        C_n = @(n) cos(n-1);
        switch model
            case 'MS_P'
                R_o = parms.R_o;
                y = R_o(end)*[1,1,1]/4;
            case {'BC_P','BCA_P'}
                L = parms.L;
                b = parms.b;
                x_s = parms.x_s;
                l_ls = parms.l_ls;
                h_s = parms.h_s;
                c = parms.c;
                arr = b:-b:-(L+2*b);
%                 y = [[x_s-l_ls*0.3,0,c]; [arr; zeros(2,numel(arr))].'];
                y = [arr; zeros(2,numel(arr))].';
%                 y(1) = -100;
            case 'S1_P'
                y = [0,0,0];
%                 y = [1,1,1]/4;
            case 'S1_P2'
                R_o = parms.R_o;
                y = R_o(end)/4*[1,1,1;
                                -1,1,1;
                                1,-1,1;
                                -1,-1,1;
                                1,1,-1;
                                -1,1,-1;
                                1,-1,-1;
                                -1,-1,-1];
            case 'Cube_P'
                a = parms.a;
                Xarr = linspace(-1,1,3)*a/4;
                Yarr = linspace(-1,1,3)*a/4;
                Zarr = linspace(-1,1,3)*a/4;
                [X,Y,Z] = ndgrid(Xarr,Yarr,Zarr);
                y = [X(:),Y(:), Z(:)];
            otherwise
                y = [0,0,0];
        end
        xms = @(v) [v(:,1)-y(1),v(:,2)-y(2),v(:,3)-y(3)];
        R = @(v) norm2(xms(v));
        analytic = @(v) analyticPulsation_(v,C_n,y,k_1);
        gAnalytic = @(v) gAnalyticPulsation_(v,C_n,y,k_1);
        dpdn = @(v,n) dAnalyticPulsation_(v,C_n,y,k_1,n);
        
        varCol{1}.farField = @(v) fAnalyticPulsation_(v,C_n,y,k_1);
        varCol{1}.analytic = analytic;
        varCol{1}.gAnalytic = gAnalytic;
        
        p_inc = @(v) zeros(size(v,1),1);
        gp_inc = @(v) zeros(size(v,1),1);
        dp_inc = @(v,n) -dpdn(v,n);
        if ~strcmp(BC,'NBC')
            error('The radial pulsation exact solution is not a solution to coupled problems')
        end
        e3Dss_options = struct('omega', omega, ...
                         'P_inc', P_inc, ...
                         'rho_f', rho_f, ...
                         'c_f', c_f);
        
    case 'planeWave'
        d_vec = -[cos(beta_s(1))*cos(alpha_s);
                  cos(beta_s(1))*sin(alpha_s);
                  sin(beta_s(1))*ones(1,length(alpha_s))];
        varCol{1}.d_vec = d_vec;
        if isSphericalShell
            R_o = parms.R_o;
            if useSolidDomain
                e3Dss_options = struct('d_vec', d_vec, ...
                                 'omega', omega, ...
                                 'R_i', R_i, ...
                                 'R_o', R_o, ...
                                 'P_inc', P_inc, ...
                                 'E', E, ...
                                 'nu', nu, ...
                                 'rho_s', rho_s, ...
                                 'rho_f', rho_f, ...
                                 'N_max', N_max,...
                                 'c_f', c_f, ...
                                 'calc_dpdx', 1,...
                                 'calc_dpdy', 1,...
                                 'calc_dpdz', 1,...
                                 'calc_u_x', 1,...
                                 'calc_u_y', 1,...
                                 'calc_u_z', 1,...
                                 'calc_du_xdx', 1,...
                                 'calc_du_xdy', 1,...
                                 'calc_du_xdz', 1,...
                                 'calc_du_ydx', 1,...
                                 'calc_du_ydy', 1,...
                                 'calc_du_ydz', 1,...
                                 'calc_du_zdx', 1,...
                                 'calc_du_zdy', 1,...
                                 'calc_du_zdz', 1,...
                                 'calc_sigma_xx', 1,...
                                 'calc_sigma_yy', 1,...
                                 'calc_sigma_zz', 1,...
                                 'calc_sigma_yz', 1,...
                                 'calc_sigma_xz', 1,...
                                 'calc_sigma_xy', 1);
            else
                e3Dss_options = struct('d_vec', d_vec, ...
                                 'omega', omega, ...
                                 'R_o', R_o, ...
                                 'P_inc', P_inc, ...
                                 'rho_f', rho_f, ...
                                 'calc_dpdx', 1,...
                                 'calc_dpdy', 1,...
                                 'calc_dpdz', 1,...
                                 'N_max', N_max,...
                                 'c_f', c_f);
            end
            varCol{1}.farField = @(v) farField_(v,e3Dss_options);
            varCol{1}.analytic = @(v) analytic_(v,e3Dss_options);
            varCol{1}.gAnalytic = @(v) gAnalytic_(v,e3Dss_options);
            if useSolidDomain
                varCol{2}.analytic = @(v) analytic_solid_(v,e3Dss_options);
            end
            if useInnerFluidDomain
                varCol{3}.analytic = @(v) analytic_fluid_i_(v,e3Dss_options);
            end
        else
            e3Dss_options = [];
        end
        if isinf(N_max)
            p_inc = @(v) P_inc*exp(1i*(v*d_vec)*k_1);
            gp_inc = @(v) 1i*elementProd(p_inc(v),d_vec.')*k_1;
            dp_inc = @(v,n) 1i*(n*d_vec)*k_1.*p_inc(v);
            dpdn = @(v,n) zeros(size(v,1),1);
        else
            p_inc = @(v) P_inc*exp(1i*(v*d_vec)*k_1);
            gp_inc = @(v) -gAnalytic_(v,e3Dss_options);
            dp_inc = @(v,n) -n*gAnalytic_(v,e3Dss_options);
            dpdn = @(v,n) zeros(size(v,1),1);
        end
end
varCol{1}.p_inc = p_inc;
varCol{1}.dp_inc = dp_inc;
varCol{1}.gp_inc = gp_inc;
if strcmp(coreMethod,'XI')
    varCol{1}.d_vec = d_vec;
end
varCol{1}.dpdn = dpdn;
varCol{1}.e3Dss_options = e3Dss_options;




