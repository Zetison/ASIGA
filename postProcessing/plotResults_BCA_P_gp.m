
for study_i = 1:numel(studies)
    study = studies(study_i);
    options = struct('xname',           'agpBEM',  ...
                     'yname',           'surfaceError', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'semilogy', ... 
                     'lineStyle',       '*-', ... 
                     'xLoopName',       'agpBEM', ...
                     'subFolderName',   '../results/BCA_P_qp', ...
                     'legendEntries',   {{'method','extraGPBEM','degree'}}, ...
                     'noXLoopPrms',     1); 
% 
    figure(2)
    printResultsToTextFiles(study,options)

%     options.yname = 'cond_number';
%     figure(3)
%     printResultsToTextFiles(study,options)
% 
%     figure(5)
%     options.legendEntries = {{'method','f','degreeElev','M'}};
%     options.noXLoopPrms = 0;
%     options.lineStyle = '-';
%     options.axisType = 'plot';
%     options.xname = 'alpha';
%     options.yname = 'TS';
%     options.xScale = 180/pi;
% %             printResultsToTextFiles(study,options)
% 
%     options.yname = 'error_p';
%     options.axisType = 'semilogy';
% 
%     figure(6)
%     printResultsToTextFiles(study,options)

%             options.yname = 'error_pAbs';
%             options.axisType = 'semilogy';
%             
%             figure(7)
%             printResultsToTextFiles(study,options)
end