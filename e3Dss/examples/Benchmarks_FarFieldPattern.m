close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models


pathToResults = '../../../results/e3Dss/';
% pathToResults = '../results';

% startMatlabPool

% Calculate Far-Field pattern (BI and Sweep)
tic
alpha_s = 240*pi/180;
beta_s = 30*pi/180;
beta_f = beta_s;
beta_f_arr = beta_s;
d_vec = -[cos(beta_s)*cos(alpha_s);
          cos(beta_s)*sin(alpha_s);
          sin(beta_s)];
models = {'S1','S3','S5','S35','S15','S13','S135'};
% models = {'S1'};

for scatteringCase = {'BI'} %,'Sweep'}
    switch scatteringCase{1}
        case 'BI'
            delta_alpha = 0.1;
            alpha_f_arr = (0:delta_alpha:360)*pi/180;
            f = [1e3, 3e3, 10e3, 30e3]; % frequencies in Hertz
%                     f = [1e3, 10e3]; % frequencies in Hertz
        case 'Sweep'
            alpha_f_arr = alpha_s;
            delta_f = 0.1;
            f = linspace(1e3,10e3,3000);
%                     f = 1e3:delta_f:30e3;
%                     f = 1e3:delta_f:5e3;
    end
    for ii = 1:length(models)
        model = models{ii};
        close all
        for SHBC = [1, 0]
            for SSBC = [1, 0]
                for ESBC = [1, 0]
                    if ~(ESBC + SHBC + SSBC > 1)
                        if strcmp(scatteringCase{1},'BI')
                            if strcmp(model,'S5') && exist([pathToResults 'S5_SHBC_F30_specialValues.mat'], 'file')
                                load([pathToResults 'S5_SHBC_F30_specialValues'])
                                alpha_f_arr = unique(sort([specialValues', alpha_f_arr]));
                            else
                                alpha_f_arr = (0:delta_alpha:360)*pi/180;
                            end
                        end
                        v = [cos(beta_f)*cos(alpha_f_arr); cos(beta_f)*sin(alpha_f_arr); sin(beta_f)*ones(size(alpha_f_arr))]'; 

                        switch model
                            case 'S1'
                                setS1Parameters
                            case 'S3'
                                setS3Parameters
                            case 'S5'
                                setS5Parameters
                            case 'S13'
                                setS13Parameters
                            case 'S15'
                                setS15Parameters
                            case 'S35'
                                setS35Parameters
                            case 'S135'
                                setS135Parameters
                            case 'tripleShell'
                                setTripleShellParameters
                        end
                        R_i = R_o - t;
                        omega = 2*pi*f; % Angular frequency

                        defineBCstring

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
                                         'calc_farField', 1);

                        % Create folders
                        if ~exist(pathToResults, 'dir')
                            mkdir(pathToResults);
                        end

                        if length(alpha_s) > 1
                            aspect = 'S';
                        else
                            aspect = num2str(round(alpha_s*180/pi, 15, 'significant'));
                        end
                        if length(beta_s) > 1
                            elevation = 'S';
                        else
                            elevation = num2str(round(beta_s*180/pi, 15, 'significant'));
                        end
                        if strcmp(scatteringCase{1}, 'Sweep')
                            frequency = 'S';
                        else
                            frequency = num2str(f/1000);
                        end
                        data = e3Dss(v, options);
                        legendArr = cell(0,1);

                        varCol.alpha_s = alpha_s;
                        varCol.beta_s = beta_s;
                        varCol.scatteringCase = scatteringCase{1};
                        if SHBC
                            colorArr = [0,70,147]/255;
                        elseif SSBC
                            colorArr = [178,0,0]/255;
                        elseif ESBC
                            colorArr = [59,124,37]/255;
                        else
                            colorArr = [149,49,157]/255;
                        end
                        if strcmp(scatteringCase{1}, 'Sweep')
                            saveName = [model '_' BC '_' scatteringCase{1} '_A'  aspect '_E' elevation '_F' frequency];
                            resultsFolderName = [pathToResults '/' saveName];
                            legendArr = {saveName};
                            TS = 20*log10(abs(data(1).p)).'; 
                            filename = [pathToResults '/' saveName];
                            varCol.saveName = saveName;
                            varCol.f_arr = f;
                            if any(~data(1).flag)
                                figure(30+ii)
                                TS_plot = TS(logical(~data(1).flag));
                                kR_01 = 2*pi*f(logical(~data(1).flag))/c_f(1)*R_o(1);
                                printResultsToFile(filename, kR_01.', TS_plot, varCol, 1, 0, 'NTNU_FFI', 'Exact solution')
                                plot(kR_01, TS_plot,'color',colorArr,'DisplayName', [model ' with ' BC]);
                                hold on
                                xlabel('$kR_{0,1}$','interpreter','latex')
                                xlim([kR_01(1), kR_01(end)])
                                legend('off');
                                legend('show','Location','northeast');
                            end
                        elseif strcmp(scatteringCase{1}, 'BI')
                            for i = 1:length(f)

                                figure(i)
                                saveName = [model '_' BC '_' scatteringCase{1} '_A'  aspect '_E' elevation '_F' num2str(f(i)/1000)];
                                filename = [pathToResults '/' saveName];
                                varCol.f = f(i);
                                varCol.saveName = saveName;
                                TS = 20*log10(abs(data(1).p(:,i))); 
                                if ~data(1).flag(i)
                                    printResultsToFile(filename, 180/pi*alpha_f_arr.', TS, varCol, 1, 0, 'NTNU_FFI', 'Exact solution')
                                    plot(alpha_f_arr*180/pi, TS,'color',colorArr,'DisplayName', [model ' with ' BC]);
                                    ylabel('TS','interpreter','latex')
                                    xlabel('$\alpha_f$','interpreter','latex')
                                    xlim([0, 360])
                                    legend('off');
                                    legend('show','Location','northeast');
                                    hold on
                                end
                            end
                        end
                    end
                end
            end
        end
        if strcmp(scatteringCase{1}, 'Sweep')
            saveName = [model '_' scatteringCase{1} '_A'  aspect '_E' elevation '_F' frequency];
            filename = [pathToResults '/' saveName];
            savefig([filename '.fig'])
        elseif strcmp(scatteringCase{1}, 'BI')
            for i = 1:length(f)
                figure(i)
                saveName = [model '_' scatteringCase{1} '_A'  aspect '_E' elevation '_F' num2str(f(i)/1000)];
                title([saveName(1:end-2) num2str(f(i)/1000)], 'interpreter', 'none')
                filename = [pathToResults '/' saveName];
                savefig([filename(1:end-2) num2str(f(i)/1000) '.fig'])
            end
        end
    end
end
toc