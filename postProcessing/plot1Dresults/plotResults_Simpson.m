close all

for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'dofs',  ...
                     'yname',           'surfaceError', ...
                     'plotResults', 	1, ... 
                     'printResults',	0, ... 
                     'axisType',        'loglog', ... 
                     'lineStyle',       '*-', ... 
                     'xLoopName',       'M', ...
                     'subFolderName',   '../results/Simpson', ...
                     'legendEntries',   {{'method','formulation','M'}}, ...
                     'noXLoopPrms',     1); 

%             figure(2)
%             printResultsToTextFiles(study,options)
% 
%             options.xname = 'tot_time';
%             options.axisType = 'loglog';
%             figure(3)
%             printResultsToTextFiles(study,options)
% 
%             options.xname = 'dofs';
%             options.yname = 'cond_number';
%             figure(4)
%             printResultsToTextFiles(study,options)

    options.noXLoopPrms = 0;
%     options.legendEntries = {'method','M','formulation'};
    options.legendEntries = {'method','M','formulation','useNeumanProj'};
    options.lineStyle = '-';
    options.xname = 'alpha';
    options.yname = 'error_pAbs';
    options.axisType = 'semilogy';
    options.xScale = 180/pi;
    options.yScale = 1/100;
    figure(5)
    printResultsToTextFiles(study,options)

    options.yname = 'abs_p';
    options.axisType = 'plot';
    options.yScale = 1;
    figure(6)
    printResultsToTextFiles(study,options)
end


figure(5)
error_simpson = importdata('../results/Simpson/imageData/Fig17_M1.csv');
loglog(180/pi*error_simpson(:,1),error_simpson(:,2),'*','DisplayName','Simpson M = 1');
error_simpson = importdata('../results/Simpson/imageData/Fig17_M2.csv');
loglog(180/pi*error_simpson(:,1),error_simpson(:,2),'*','DisplayName','Simpson M = 2');
error_simpson = importdata('../results/Simpson/imageData/Fig17_M3.csv');
loglog(180/pi*error_simpson(:,1),error_simpson(:,2),'*','DisplayName','Simpson M = 3');



