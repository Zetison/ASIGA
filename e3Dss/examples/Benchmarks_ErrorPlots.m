close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

pathToResults = '../../plotData/e3Dss/';
% pathToResults = '../results';
% mpstartup
startMatlabPool
%% Calculate errors
for useSymbolicPrecision = [0,1]
    if useSymbolicPrecision
        prec = 'mp';
%         prec = 'sym';
    else
        prec = 'double';
    end
    models = {'S1','S3','S5','S35','S15','S13','S135'};
%     models = {'S1'};
    counter = 1;
    for i_model = 1:length(models)
        for ESBC = [0, 1]
            for SHBC = [0, 1]
                for SSBC = [0, 1]
                    if ~(ESBC + SHBC + SSBC > 1)
                        tasks(counter).model = models{i_model};
                        tasks(counter).ESBC = ESBC;
                        tasks(counter).SHBC = SHBC;
                        tasks(counter).SSBC = SSBC;
                        counter = counter + 1;
                    end
                end
            end
        end
    end
%     for i = 1:length(models)
    parfor i = 1:length(tasks)
        switch prec
            case 'single'
                Eps = 1e-7;
                PI = pi;
                O = 0;
            case 'double'
                Eps = eps;
                PI = pi;
                O = 0;
            case 'sym'
                Eps = 1e-50;
                digits(2000);
                PI = vpa('pi');
                O = vpa('0');
            case 'mp'
                Eps = 1e-50;
                mp.Digits(2000);
                PI = mp('pi');
                O = mp('0');
        end
        model = tasks(i).model;
        ESBC = tasks(i).ESBC;
        SHBC = tasks(i).SHBC;
        SSBC = tasks(i).SSBC;
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
        if strcmp(prec,'sym')
            P_inc = vpa(P_inc);
            rho_s = vpa(rho_s);
            rho_f = vpa(rho_f);
            c_f = vpa(c_f);
            t = vpa(t);
            E = vpa(E);
            nu = vpa(nu);
            R_o = vpa(R_o); % Density of solid
        elseif strcmp(prec,'mp')
            P_inc = mp(P_inc);
            rho_s = mp(rho_s);
            rho_f = mp(rho_f);
            c_f = mp(c_f);
            t = mp(t);
            E = mp(E);
            nu = mp(nu);
            R_o = mp(R_o); % Density of solid

        end

        R_i = R_o - t; % Inner radius of shell
        alpha_s = 240*PI/180;
        beta_s = 30*PI/180;
        d_vec = zeros(3,1,class(PI));
        d_vec(1) = -cos(beta_s)*cos(alpha_s);
        d_vec(2) = -cos(beta_s)*sin(alpha_s);
        d_vec(3) = -sin(beta_s);

%                             defineBCstring
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
                r = linspaceHP(R_o(1), 2*R_o(1), npts_r);
            elseif mod(j,2)
                if m == M+1
                    r = linspaceHP(O, R_i(m-1), npts_r);
                else
                    r = linspaceHP(R_o(m), R_i(m-1), npts_r);
                end
            else
                if ESBC && m == M
                    r = linspaceHP(O, R_o(m), npts_r);    
                else
                    r = linspaceHP(R_i(m), R_o(m), npts_r);     
                end
                m = m + 1;
            end
            theta = linspaceHP(O,PI,npts_theta);
            phi = linspaceHP(O,2*PI,npts_phi);
            pts = zeros(length(r)*length(theta)*length(phi),3,class(PI));

            counter = 1;
            for ii = 1:length(r)
                for jj = 1:length(theta)
                    for ll = 1:length(phi)
                        pts(counter,:) = r(ii)*[sin(theta(jj))*cos(phi(ll)), sin(theta(jj))*sin(phi(ll)), cos(theta(jj))];
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

        C = (R_o(1)./c_f(1))^(3/2)/Upsilon^(1/2);
        if 1
            nFreqs = 2; %
%             f = 10.^linspaceHP(-log10(1e3*C),-log10(5e2*C),nFreqs);
            f = 10.^linspaceHP(-log10(1e3*C),log10(4e2/C),nFreqs);
        else
            nFreqs = 100; %
            f = 10.^linspaceHP(-log10(1e3*C),log10(4e2/C),nFreqs);
        end

        omega = 2*PI*f;
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
                         'calc_u_x', 1, ...
                         'calc_u_y', 1, ...
                         'calc_u_z', 1, ...
                         'calc_dpdx', 1, ...
                         'calc_dpdy', 1, ...
                         'calc_dpdz', 1, ...
                         'calc_p_laplace', 1, ...
                         'calc_sigma_xx', 1, ...
                         'calc_sigma_yy', 1, ...
                         'calc_sigma_zz', 1, ...
                         'calc_sigma_yz', 1, ...
                         'calc_sigma_xz', 1, ...
                         'calc_sigma_xy', 1, ...
                         'calc_sigma_rr', 1, ...
                         'calc_navier', 1, ...
                         'calc_errorsNavier', 1, ...
                         'calc_errorsDisplacementCondition', 1, ...
                         'calc_errorsPressureCondition', 1, ...
                         'calc_errorsHelmholtz', 1, ...
                         'useSymbolicPrecision', useSymbolicPrecision, ...
                         'prec', prec, ...
                         'Eps', Eps);

        data = e3Dss(vv, options);
        err_navier1 = zeros(M,nFreqs,class(PI));
        err_navier2 = zeros(M,nFreqs,class(PI));
        err_helmholtz = zeros(M+1,nFreqs,class(PI));
        err_pc = zeros(M,nFreqs,class(PI));
        err_dc = zeros(M,nFreqs,class(PI));
        for m = 1:M+1
            if m ~= M+1
                err_dc(m,:) = data(m).err_dc;
                if ~(SHBC && m == M)
                    err_pc(m,:) = data(m).err_pc;
                    err_navier1(m,:) = data(m).err_navier1;
                    err_navier2(m,:) = data(m).err_navier2;
                end
            end
            if ~(m == M+1 && (SHBC || SSBC || ESBC))
                err_helmholtz(m,:) = data(m).err_helmholtz;
            end
        end
        err_navier1 = max(err_navier1,[],1);
        err_navier2 = max(err_navier2,[],1);
        err_helmholtz = max(err_helmholtz,[],1);
        err_pc = max(err_pc,[],1);
        err_dc = max(err_dc,[],1);

        sc = f*C;
        figure
        if M == 1 && SHBC
            loglog(sc, err_helmholtz,'color',[0,70,147]/255,'DisplayName','Helmholtz')
            hold on
            loglog(sc, err_dc,'color',[178,0,0]/255,'DisplayName','Displacenemnt condition')
            legendArr = {'Helmholtz', 'DisplacementCond'};
            err = [err_helmholtz; err_dc];
        else
            loglog(sc, err_helmholtz,'color',[0,70,147]/255,'DisplayName','Helmholtz')
            hold on
            loglog(sc, err_navier1,'color',[178,0,0]/255,'DisplayName', 'Navier - $1^{\mathrm{st}}$ component')
            loglog(sc, err_navier2,'color',[59,124,37]/255,'DisplayName','Navier - $2^{\mathrm{nd}}$ component')
            loglog(sc, err_dc,'color',[149,49,157]/255,'DisplayName','Displacenemnt condition')
            loglog(sc, err_pc,'color',[247, 158,30]/255,'DisplayName','Pressure condition')
            legendArr = {'Helmholtz', 'Navier1', 'Navier2', 'DisplacementCond', 'PressureCond'};
            err = [err_helmholtz; err_navier1; err_navier2; err_dc; err_pc];
        end
        leg1 = legend('show','Location','northwest');
        set(leg1,'Interpreter','latex');
        filename = [pathToResults 'errors_' model '_' BC '_Symbolic' num2str(useSymbolicPrecision)];

        printResultsToFile(filename, double(sc.'), double(err.'), [], 0, 1, [], [], {'Cf'},legendArr)
        xlabel('$C f$','interpreter','latex')
        ylabel('Relative residual error')
        title(['Errors for model ' model '_' BC], 'interpreter', 'none')
        hold off
        if ~useSymbolicPrecision
            ylim([0.1*eps 1e2])
        end
        xlim([double(sc(1)), double(sc(end))])
        drawnow
        savefig([filename '.fig'])
        fprintf('Finished a case in %f seconds!\n\n', toc)
    end
end

