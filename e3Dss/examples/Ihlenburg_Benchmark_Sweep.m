close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

pathToResults = '../../../results/e3Dss/';
% pathToResults = '../results';

startMatlabPool

ESBC = 0;
SSBC = 0;
k_arr = linspace(0.001, 2, 3000).';
k_arr = [];
specialValues = [0.250621182794531 %
               0.320300579445871 %
               0.370671479527136 %
               0.412992731010227 %
               0.454191270410376 %
               0.499088778976889 %
               0.551286412239756 %
               0.613370456080303 %
               0.687008309336546 %
               0.773084257718564 %
               0.871890313908958 %
               0.983323027396819 %
               1.107045032953710 %
               1.242597693362253 %
               1.389470517759271 %
               1.547139666101034 %
               1.715087015757087 %
               1.892808062465205 %
               ];
k_arr = unique(sort([k_arr;specialValues]));
for SHBC = 0 %[0, 1]
    for modelCell = {'IL'} %{'IL', 'S5', 'S35', 'S135'}
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
        parfor i = 1:length(k_arr)
            k = k_arr(i);
            omega = k*c_f(1);
            R_a = 2*R_o(1);
            
            alpha_s = pi;
            beta_s = 0;
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
            vtfFileName = [folderName BC '_' num2str(i)];

            extraPts = 15;

            createParaviewFiles_exact3(extraPts, vtfFileName, options)
        end
    end
end

% %         solid = getSphericalShellData(5, 4.992,'Zaxis');
% %         solid = getSphericalShellData(3, 2.98,'Zaxis');
%         solid = getSphericalShellData(1, 0.95,'Zaxis');
%         [nodes, ~, visElements] = buildVisualization3dMesh_new3(solid.knots{1}, solid.knots{2}, solid.knots{3}, 200, 200, 0, solid);
%         VTKoptions = struct('name','otherFunctions/e3Dss/results/sources/S1_solid', 'celltype', 'VTK_HEXAHEDRON'); 
%         VTKdata.nodes = nodes;
%         VTKdata.visElements = visElements;
%         VTKdata.omega = 1;
%         makeVTKfile(VTKdata, VTKoptions);