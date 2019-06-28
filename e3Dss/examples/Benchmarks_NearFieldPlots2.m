close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

pathToResults = '../../../results/e3Dss/';
% pathToResults = '../results';

startMatlabPool
N_max = 2;
ESBC = 0;
SSBC = 0;
for SHBC = 1 %[0, 1]
    for modelCell = {'S1'} %{'IL', 'S5', 'S35', 'S135'}
        model = modelCell{1};
        switch model
            case 'S1'
                setS1Parameters
            case 'S5'
                setS5Parameters
            case 'S35'
                setS35Parameters
            case 'S135'
                setS135Parameters
            case 'IL'
                setIhlenburgParameters
        end
        R_i = R_o - t;
        defineBCstring
%         k_arr = [1.892808062465205, 1.8929];
        plotRealPart = 1;
        for f = 10000
            omega = 2*pi*f;
            % rho_f(end) = 0;
            R_a = 2; %R_o(1);
            
            alpha_s = 240*pi/180;
            beta_s = 30*pi/180;
            beta_f = beta_s;
            beta_f_arr = beta_s;
            d_vec = -[cos(beta_s)*cos(alpha_s);
                      cos(beta_s)*sin(alpha_s);
                      sin(beta_s)]; 
%             d_vec = [0, 0, 1]'; 
                  
            usePointChargeWave = false;
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
                             'SHBC', SHBC, ...
                             'ESBC', ESBC, ...
                             'SSBC', SSBC, ...
                             'N_max', N_max, ...
                             'computeForSolidDomain', 0, ...
                             'plotTimeOscillation', 0, ...
                             'plotInTimeDomain', 0, ...
                             'usePointChargeWave', usePointChargeWave, ...
                             'usePlaneWave', ~usePointChargeWave, ...
                             'R_a', R_a);

            folderName = [pathToResults 'paraviewResults/' model '/'];
%             folderName = ['results/paraviewResults/' model '/'];
            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end
            vtfFileName = [folderName BC '_f_' num2str(f) '_N_' num2str(N_max)];

            extraPts = 160;

            createParaviewFiles_exact3(extraPts, vtfFileName, options)
        end
    end
end
% options.name = vtfFileName;
% solid = getSphericalShellData(1, 0.95,'Zaxis');
% artificialBoundary = extractOuterSurface(solid);
% plotModelInParaview(artificialBoundary{1}, 200, 200, 0, options, 0, 'scatterer')
% solid = getSphericalShellData(1+2*pi/(32-pi), 0.95,'Zaxis');
% artificialBoundary = extractOuterSurface(solid);
% plotModelInParaview(artificialBoundary{1}, 200, 200, 0, options, 0, 'artificialBoundary')

        