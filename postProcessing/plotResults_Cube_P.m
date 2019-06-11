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
                     'subFolderName',   '../results/Cube_P', ...
                     'legendEntries',   {{'method','formulation','degree','useNeumanProj','colBEM_C0'}}, ...
                     'noXLoopPrms',     1); 

    figure(2)
    printResultsToTextFiles(study,options)
end