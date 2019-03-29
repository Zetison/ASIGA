close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

pathToResults = '../../plotData/e3Dss/';
% pathToResults = '../results';

startMatlabPool

%% Calculate errors
useSymbolicPrecision = 0;
if useSymbolicPrecision
    Eps = 1e-40;
else
    Eps = eps;
end

models = {'S1','S3','S5','S35','S15','S13','S135'};

poolobj = gcp('nocreate');
noCoresToUse = feature('numCores');
if noCoresToUse > 7
    noCoresToUse = 7;
end
if isempty(poolobj)
    parpool(noCoresToUse, 'IdleTimeout', Inf)
end
tic
noTestCases = 1;
for i = 1:length(models)
% parfor i = 1:length(models)
    model = models{i};
    for ESBC = [0, 1]
        for SHBC = [0, 1]
            for SSBC = [0, 1]
                if ~(ESBC + SHBC + SSBC > 1)
                    tic
                    switch model
                        case 'S1'
%                                     setS1Parameters
                            P_inc = 1; % Amplitude of incident wave
                            rho_s = 7850; % Density of solid
                            rho_f = [1000, 1.2]; % Density of fluids
                            c_f = [1500, 340];  % Speed of sound in fluid domains
                            t = 0.05; % The thickness of the sphere
                            E = 210e9; % Youngs modulus of elastic material
                            nu = 0.3; % Poisson ratio of elastic material

                            R_o = 1; % Distance to far field point

                        case 'S3'
%                                     setS3Parameters
                            P_inc = 1; % Amplitude of incident wave
                            rho_s = 7850; % Density of solid
                            rho_f = [1000, 1.2]; % Density of fluids
                            c_f = [1500, 340];  % Speed of sound in fluid domains
                            t = 0.02; % The thickness of the sphere
                            E = 210e9; % Youngs modulus of elastic material
                            nu = 0.3; % Poisson ratio of elastic material

                            R_o = 3; % Distance to far field point
                        case 'S5'
%                                     setS5Parameters
                            P_inc = 1; % Amplitude of incident wave
                            rho_s = 7850; % Density of solid
                            rho_f = [1000, 1000]; % Density of fluids
                            c_f = [1500, 1500];  % Speed of sound in fluid domains
                            t = 0.008; % The thickness of the sphere
                            E = 210e9; % Youngs modulus of elastic material
                            nu = 0.3; % Poisson ratio of elastic material

                            R_o = 5; % Distance to far field point
                        case 'S13'
%                                     setS13Parameters
                            P_inc = 1; % Amplitude of incident wave
                            rho_s = [7850, 7850]; % Density of solid
                            rho_f = [1000, 1.2, 1.2]; % Density of fluids
                            c_f = [1500, 340, 340];  % Speed of sound in fluid domains
                            t = [0.02 0.05]; % The thickness of the sphere
                            E = [210e9, 210e9]; % Youngs modulus of elastic material
                            nu = [0.3, 0.3]; % Poisson ratio of elastic material

                            R_o = [3, 1];
                        case 'S15'
%                                     setS15Parameters
                            P_inc = 1; % Amplitude of incident wave
                            rho_s = [7850, 7850]; % Density of solid
                            rho_f = [1000, 1000, 1.2]; % Density of fluids
                            c_f = [1500, 1500, 340];  % Speed of sound in fluid domains
                            t = [0.008, 0.05]; % The thickness of the sphere
                            E = [210e9, 210e9]; % Youngs modulus of elastic material
                            nu = [0.3, 0.3]; % Poisson ratio of elastic material

                            R_o = [5, 1];

                        case 'S35'
%                                     setS35Parameters
                            P_inc = 1; % Amplitude of incident wave
                            rho_s = [7850, 7850]; % Density of solid
                            rho_f = [1000, 1000, 1.2]; % Density of fluids
                            c_f = [1500, 1500, 340];  % Speed of sound in fluid domains
                            t = [0.008, 0.02]; % The thickness of the sphere
                            E = [210e9, 210e9]; % Youngs modulus of elastic material
                            nu = [0.3, 0.3]; % Poisson ratio of elastic material

                            R_o = [5, 3]; 

                        case 'S135'
%                                     setS135Parameters
                            P_inc = 1; % Amplitude of incident wave
                            rho_s = [7850, 7850, 7850]; % Density of solid
                            rho_f = [1000, 1000, 1.2, 1.2]; % Density of fluids
                            c_f = [1500, 1500, 340, 340];  % Speed of sound in fluid domains
                            t = [0.008, 0.02, 0.05]; % The thickness of the sphere
                            E = [210e9, 210e9, 210e9]; % Youngs modulus of elastic material
                            nu = [0.3, 0.3, 0.3]; % Poisson ratio of elastic material

                            R_o = [5, 1, 0.5];

                    end
                    R_i = R_o - t; % Inner radius of shell
                    alpha_s = 240*pi/180;
                    beta_s = 30*pi/180;

                    if useSymbolicPrecision
                        P_inc = vpa(P_inc);
                        rho_s = vpa(rho_s);
                        rho_f = vpa(rho_f);
                        c_f = vpa(c_f);
                        t = vpa(t);
                        E = vpa(E);
                        nu = vpa(nu);
                        R_o = vpa(R_o); % Density of solid
                        R_i = R_o - t;
                        alpha_s = 240*vpa(pi)/180;
                        beta_s = 30*vpa(pi)/180;
                    end
                    d_vec = -[cos(beta_s)*cos(alpha_s);
                              cos(beta_s)*sin(alpha_s);
                              sin(beta_s)];

%                             defineBCstring
                    M = length(R_o);
                    if SHBC
                        E = E(1:end-1);
                        rho_s = rho_s(1:end-1);
                        rho_f = rho_f(1:end-1);
                        c_f = c_f(1:end-1);
                        nu = nu(1:end-1);
                        R_i = R_i(1:end-1);
                        BC = 'SHBC';
                    elseif ESBC
                        rho_f = rho_f(1:end-1);
                        c_f = c_f(1:end-1);
                        R_i = R_i(1:end-1);
                        BC = 'ESBC';
                    elseif SSBC
                        rho_f = rho_f(1:end-1);
                        c_f = c_f(1:end-1);
                        BC = 'SSBC';
                    else
                        BC = 'NNBC';
                    end

                    M = length(R_o);
                    noDomains = 2*M+1-ESBC-2*SHBC-SSBC;
                    vv = cell(noDomains,1);
                    m = 1;
                    npts_r = 4;
                    npts_theta = 4;
                    npts_phi = 4;
                    for j = 1:noDomains
                        if j == 1
                            r = linspace(R_o(1), 2*R_o(1), npts_r);
                        elseif mod(j,2)
                            if m == M+1
                                r = linspace(0, R_i(m-1), npts_r-1);
                            else
                                r = linspace(R_o(m), R_i(m-1), npts_r);
                            end
                        else
                            if ESBC && m == M
                                r = linspace(0, R_o(m), npts_r-1);    
                            else
                                r = linspace(R_i(m), R_o(m), npts_r);     
                            end
                            m = m + 1;
                        end
                        theta = linspace(0,pi,npts_theta);
                        phi = linspace(0,2*pi,npts_phi);
                        pts = zeros(length(r)*length(theta)*length(phi),3);
                        if useSymbolicPrecision
                            r = vpa(r);
                            theta = vpa(theta);
                            phi = vpa(phi);
                            pts = vpa(pts);
                        end
                        counter = 1;
                        for ii = 1:length(r)
                            for jj = 1:length(theta)
                                for ll = 1:length(phi)
                                    pts(counter,:) = r(ii)*[sin(theta(jj)).*cos(phi(ll)), sin(theta(jj)).*sin(phi(ll)), cos(theta(jj))];
                                    counter = counter + 1;
                                end
                            end
                        end
                        [~, I, ~] = uniquetol(double(pts),10*eps,'ByRows',true, 'DataScale',max(max(abs(double(pts)))));
                        vv{j} = pts(I,:);
                    end
                    K = E./(3*(1-2*nu));
                    G = E./(2*(1+nu));
                    c_s_1 = sqrt((3*K+4*G)./(3*rho_s));
                    c_s_2 = sqrt(G./rho_s);

                    switch BC
                        case {'ESBC', 'SHBC'}
                            Upsilon = min([R_i./c_s_1(1:end-1), R_i./c_s_2(1:end-1), R_o./c_f]);
                        case 'SSBC'
                            Upsilon = min([R_i./c_s_1, R_i./c_s_2, R_o./c_f]);
                        case 'NNBC'
                            Upsilon = min([R_i./c_s_1, R_i./c_s_2, R_o./c_f(1:end-1)]);
                    end
                    
                    C = (R_o(1)./c_f(1))^(3/2) ...
                            /Upsilon^(1/2);
                    nFreqs = 100;
                    f = 10.^linspace(log10(1e-3/C),log10(4e2/C),nFreqs);
                    if useSymbolicPrecision
                        f = vpa(f);
                    end
                    omega = 2*pi*f;
                    options = struct('d_vec', d_vec, ...
                                     'omega', omega, ...
                                     'R_i', R_i, ...
                                     'R_o', R_o, ...
                                     'P_inc', P_inc, ...
                                     'E', E, ...
                                     'nu', nu, ...
                                     'rho_s', rho_s, ...
                                     'rho_f', rho_f, ...
                                     'c_f', c_f, ...
                                     'calc_du_xdx', 1, ...
                                     'calc_du_xdy', 1, ...
                                     'calc_du_xdz', 1, ...
                                     'calc_du_ydx', 1, ...
                                     'calc_du_ydy', 1, ...
                                     'calc_du_ydz', 1, ...
                                     'calc_du_zdx', 1, ...
                                     'calc_du_zdy', 1, ...
                                     'calc_du_zdz', 1, ...
                                     'calc_sigma_xx', 1, ...
                                     'calc_sigma_yy', 1, ...
                                     'calc_sigma_zz', 1, ...
                                     'calc_sigma_yz', 1, ...
                                     'calc_sigma_xz', 1, ...
                                     'calc_sigma_xy', 1, ...
                                     'useSymbolicPrecision', useSymbolicPrecision, ...
                                     'Eps', Eps);

                    options.useSymbolicPrecision = useSymbolicPrecision;
                    data = e3Dss(vv, options);
                    
                    err_stress_xx = zeros(M,nFreqs);
                    err_stress_yy = zeros(M,nFreqs);
                    err_stress_zz = zeros(M,nFreqs);
                    err_stress_yz = zeros(M,nFreqs);
                    err_stress_xz = zeros(M,nFreqs);
                    err_stress_xy = zeros(M,nFreqs);
                    for m = 1:M+1
                        if m ~= M+1
                            if ~(SHBC && m == M)                          
                                strain = cell(6,1);
                                for ii = 1:6
                                    strain{ii} = zeros(size(vv{2*m},1),nFreqs);
                                    if useSymbolicPrecision
                                        strain{ii} = vpa(strain{ii});
                                    end
                                end                         
                                stress = strain;             
                                du_X = cell(3,3);
                                du_X{1,1} = data(m).du_xdx;
                                du_X{1,2} = data(m).du_xdy;
                                du_X{1,3} = data(m).du_xdz;
                                du_X{2,1} = data(m).du_ydx;
                                du_X{2,2} = data(m).du_ydy;
                                du_X{2,3} = data(m).du_ydz;
                                du_X{3,1} = data(m).du_zdx;
                                du_X{3,2} = data(m).du_zdy;
                                du_X{3,3} = data(m).du_zdz;
                                vgtinv = [1 6 5;
                                          6 2 4;
                                          5 4 3];
                                for ii = 1:3
                                    for jj = ii:3
                                        strain{vgtinv(ii,jj)} = 0.5*(du_X{ii,jj}+du_X{jj,ii});
                                    end
                                end
                                for ii = 1:3
                                    for jj = 1:3
                                        if ii == jj
                                            stress{ii} = stress{ii} + (K(m)+4*G(m)/3)*strain{jj};
                                        else
                                            stress{ii} = stress{ii} + (K(m)-2*G(m)/3)*strain{jj};
                                        end
                                    end
                                end
                                for ii = 4:6
                                    stress{ii} = 2*G(m)*strain{ii};
                                end
                                err_stress_xx(m,:) = max(abs(stress{1} - data(m).sigma_xx),[],1)./max(abs(data(m).sigma_xx),[],1);
                                err_stress_yy(m,:) = max(abs(stress{2} - data(m).sigma_yy),[],1)./max(abs(data(m).sigma_yy),[],1);
                                err_stress_zz(m,:) = max(abs(stress{3} - data(m).sigma_zz),[],1)./max(abs(data(m).sigma_zz),[],1);
                                err_stress_yz(m,:) = max(abs(stress{4} - data(m).sigma_yz),[],1)./max(abs(data(m).sigma_yz),[],1);
                                err_stress_xz(m,:) = max(abs(stress{5} - data(m).sigma_xz),[],1)./max(abs(data(m).sigma_xz),[],1);
                                err_stress_xy(m,:) = max(abs(stress{6} - data(m).sigma_xy),[],1)./max(abs(data(m).sigma_xy),[],1);
                            end
                        end
                    end
                    err_stress_xx = max(err_stress_xx,[],1);
                    err_stress_yy = max(err_stress_yy,[],1);
                    err_stress_zz = max(err_stress_zz,[],1);
                    err_stress_yz = max(err_stress_yz,[],1);
                    err_stress_xz = max(err_stress_xz,[],1);
                    err_stress_xy = max(err_stress_xy,[],1);

                    sc = f*C;
                    figure(noTestCases)
                    if ~(M == 1 && SHBC)
                        loglog(sc, err_stress_xx,'color',[0,70,147]/255,'DisplayName', '$\sigma_{xx}$')
                        hold on
                        loglog(sc, err_stress_yy,'color',[178,0,0]/255,'DisplayName', '$\sigma_{yy}$')
                        loglog(sc, err_stress_zz,'color',[59,124,37]/255,'DisplayName', '$\sigma_{zz}$')
                        loglog(sc, err_stress_yz,'color',[149,49,157]/255,'DisplayName', '$\sigma_{yz}$')
                        loglog(sc, err_stress_xz,'color',[247, 158,30]/255,'DisplayName', '$\sigma_{xz}$')
                        loglog(sc, err_stress_xy,'color',[0,172,239]/255,'DisplayName', '$\sigma_{xy}$')
                        err = [err_stress_xx; err_stress_yy; err_stress_zz; err_stress_yz; err_stress_xz; err_stress_xy];
                    else
                        continue
                    end
                    leg1 = legend('show','Location','northwest');
                    set(leg1,'Interpreter','latex');
                    filename = [pathToResults 'sigmaErrors_' model '_' BC '_Symbolic' num2str(useSymbolicPrecision)];

%                     printResultsToFile(filename, double(sc.'), double(err.'), [], 0, 1, [], [], {'Cf'},legendArr)
                    xlabel('$C f$','interpreter','latex')
                    ylabel('Relative error')
                    title(['Errors for model ' model '_' BC], 'interpreter', 'none')
                    hold off
                    if ~useSymbolicPrecision
                        ylim([0.1*eps 1e2])
                    end
                    xlim([double(sc(1)), double(sc(end))])
                    drawnow
%                     savefig([filename '.fig'])
                    noTestCases = noTestCases + 1;
                    fprintf('Finished a case in %f seconds!\n\n', toc)
                end
            end
        end
    end
end

