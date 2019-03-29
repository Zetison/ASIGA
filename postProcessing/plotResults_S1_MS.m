
for study_i = 1:numel(studies)
    study = studies(study_i);
    options = struct('xname',           'nepw',  ...
                     'yname',           'surfaceError', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'loglog', ... 
                     'lineStyle',       '*-', ... 
                     'xLoopName',       'M', ...
                     'subFolderName',   'results/_studies/S1_MS', ...
                     'legendEntries',   {{'method','formulation','N'}}, ...
                     'noXLoopPrms',     1); 
% 
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
% 
    options.noXLoopPrms = 0;
    options.legendEntries = {'method','M','parm'};
    options.lineStyle = '-';
    options.xname = 'alpha';
    options.yname = 'error_p';
    options.axisType = 'semilogy';
    options.xScale = 180/pi;
    figure(6)
    printResultsToTextFiles(study,options)

    options.yname = 'error_pAbs';
    figure(7)
    printResultsToTextFiles(study,options)

    options.yname = 'TS';
    options.axisType = 'plot';
    options.yScale = 1;
    figure(5)
    printResultsToTextFiles(study,options)
end