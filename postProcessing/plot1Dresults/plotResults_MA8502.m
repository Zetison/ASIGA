
for study_i = 1:numel(studies)
    study = studies(study_i);
    options = struct('xname',           'dofsAlg',  ... % dofsAlg
                     'yname',           'energyError', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'semilogy', ... 
                     'lineStyle',       '*-', ... 
                     'xLoopName',       'degreeElev', ...
                     'subFolderName',   'results/MA8502', ...
                     'legendEntries',   {{'coreMethod','f','N'}}, ...
                     'noXLoopPrms',     1); 
% 
    figure(2)
    printResultsToTextFiles(study,options)

    figure(3)
    options.axisType = 'loglog';
    options.xname = 'dofs';
    options.yname = 'cond_number';
    printResultsToTextFiles(study,options)

    figure(4)
    options.yname = 'energyError';
    options.xname = 'tot_time';
    printResultsToTextFiles(study,options)

    figure(5)
    options.xname = 'timeSolveSystem';
    printResultsToTextFiles(study,options)

    figure(6)
    options.xname = 'timeBuildSystem';
    printResultsToTextFiles(study,options)
end