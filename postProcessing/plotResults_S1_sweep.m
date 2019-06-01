
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'f',  ...
                     'yname',           'TS', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'plot', ... 
                     'lineStyle',       '-', ... 
                     'xLoopName',       'f', ...
                     'subFolderName',   '../results/_studies/S1_sweep/', ...
                     'legendEntries',   {{'method','M','parm'}}, ...
                     'noXLoopPrms',     1); 

    figure(2)
    printResultsToTextFiles(study,options)

    options.yname = 'error_p';
    options.axisType = 'semilogy';

    figure(4)
    printResultsToTextFiles(study,options)

    options.yname = 'error_pAbs';
    options.axisType = 'semilogy';

    figure(5)
    printResultsToTextFiles(study,options)

    options.yname = 'surfaceError';

    figure(6)
    printResultsToTextFiles(study,options)
end