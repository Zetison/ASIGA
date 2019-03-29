
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'dofs',  ...
                     'yname',           'surfaceError', ...
                     'plotResults', 	1, ... 
                     'printResults',	0, ... 
                     'axisType',        'loglog', ... 
                     'lineStyle',       '*-', ... 
                     'xLoopName',       'M', ...
                     'noXLoopPrms',     1); 

    figure(2)
    printResultsToTextFiles(study,options)
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
end