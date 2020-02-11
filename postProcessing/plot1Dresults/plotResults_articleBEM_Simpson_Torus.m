close all
for study_i = 1:numel(studies)
    study = studies(study_i);
    options = struct('xname',           'nepw',  ...
                     'yname',           'surfaceError', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'loglog', ... 
                     'lineStyle',       '*-', ... 
                     'xLoopName',       'M', ...
                     'yScale',          1/100, ... 
                     'subFolderName',   '../results/articleBEM_Simpson_Torus', ...
                     'legendEntries',   {{'method','coreMethod','formulation','degree','extraGP','extraGPBEM','agpBEM','useNeumanProj'}}, ...
                     'noXLoopPrms',     1); 

    figure(2)
    printResultsToTextFiles(study,options)

%     figure(3)
    options.yname = 'energyError';
    printResultsToTextFiles(study,options)

    options.yname = 'surfaceError';
    options.xname = 'h_max';
    options.yScale = 1/100;
    
%     figure(4)
%     printResultsToTextFiles(study,options)
end

% figure(2)
% error_simpson = importdata('../results/Simpson_Torus/imageData/Fig24_p2.csv');
% loglog(error_simpson(:,1),error_simpson(:,2),'*-','DisplayName','Simpson p = 2');
% error_simpson = importdata('../results/Simpson_Torus/imageData/Fig24_p3.csv');
% loglog(error_simpson(:,1),error_simpson(:,2),'*-','DisplayName','Simpson p = 3');

% figure(4)
% error_simpson = importdata('../results/articleBEM_Simpson_Torus/imageData/Fig23_p2_mod.csv');
% loglog(error_simpson(:,1),error_simpson(:,2),'*-','DisplayName','Simpson p = 2');
% error_simpson = importdata('../results/articleBEM_Simpson_Torus/imageData/Fig23_p3_mod.csv');
% loglog(error_simpson(:,1),error_simpson(:,2),'*-','DisplayName','Simpson p = 3');
% error_simpson = importdata('../results/articleBEM_Simpson_Torus/imageData/Fig23_p4_mod.csv');
% loglog(error_simpson(:,1),error_simpson(:,2),'*-','DisplayName','Simpson p = 4');
% error_simpson = importdata('../results/articleBEM_Simpson_Torus/imageData/Fig23_p5_mod.csv');
% loglog(error_simpson(:,1),error_simpson(:,2),'*-','DisplayName','Simpson p = 5');
