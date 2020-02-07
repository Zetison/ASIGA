close all
for study_i = 1:numel(studies)
    study = studies(study_i);
    options = struct('xname',           'M',  ...
                     'yname',           'surfaceError', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'semilogy', ... 
                     'lineStyle',       '*-', ... 
                     'xLoopName',       'M', ...
                     'legendEntries',   {{'method'}}, ...
                     'noXLoopPrms',     1); 

    figure(2)
    printResultsToTextFiles(study,options)
    
    figure(3)
    options.xname = 'alpha';
    options.yname = 'TS';
    options.lineStyle = '-';
    options.xScale = 180/pi;
    options.noXLoopPrms = 0;
    options.axisType = 'plot';
    options.legendEntries = {'method','M'};
    printResultsToTextFiles(study,options)
    
    figure(4)
    options.axisType = 'polar';
    printResultsToTextFiles(study,options)
end