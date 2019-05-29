
close all
for study_i = 1:numel(studies)
    study = studies(study_i);
    options = struct('xname',           'agpBEM',  ...
                     'yname',           'surfaceError', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'semilogy', ... 
                     'lineStyle',       '*-', ... 
                     'xLoopName',       'agpBEM', ...
                     'subFolderName',   '../results/S1_BEM_qp', ...
                     'legendEntries',   {{'extraGP','extraGPBEM','colBEM_C0','formulation','quadMethodBEM','degree','M'}}, ...
                     'noXLoopPrms',     1); 

    figure(2)
    printResultsToTextFiles(study,options)
% 
%     options.xname = 'nepw';
%     options.yScale = 1;
%     figure(3)
%     printResultsToTextFiles(study,options)

%     options.xname = 'h_max';
%     options.yScale = 1/100;
% %     options.xScale = 0.2874/3.7417;
%     figure(4)
%     printResultsToTextFiles(study,options)
% 
%             options.xname = 'dofs';
%             options.yname = 'cond_number';
%             figure(4)
%             printResultsToTextFiles(study,options)
end
