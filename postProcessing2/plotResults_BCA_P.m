
for study_i = 1:numel(studies)
    study = studies(study_i);
    options = struct('xname',           'dofs',  ...
                     'yname',           'surfaceError', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'loglog', ... 
                     'lineStyle',       '*-', ... 
                     'xLoopName',       'M', ...
                     'subFolderName',   '../results/BCA_P', ...
                     'legendEntries',   {{'method','formulation','M','degree'}}, ...
                     'noXLoopPrms',     1); 

    figure(2)
    printResultsToTextFiles(study,options)
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
    options.legendEntries = {'method','M','formulation'};
    options.lineStyle = '-';
    options.xname = 'alpha';
    options.yname = 'error_pAbs';
    options.axisType = 'semilogy';
    options.xScale = 180/pi;
%     options.yScale = 1/100;
    figure(5)
    printResultsToTextFiles(study,options)

    options.yname = 'TS';
    options.axisType = 'plot';
    options.yScale = 1;
    figure(6)
    printResultsToTextFiles(study,options)
end