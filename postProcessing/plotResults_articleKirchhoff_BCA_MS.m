close all
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'alpha',  ...
                     'yname',           'TS', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'plot', ... 
                     'lineStyle',       '-', ... 
                     'subFolderName',   '../results/articleKirchhoff_BCA_MS', ...
                     'legendEntries',   {{'method','M','coreMethod'}}, ...
                     'noXLoopPrms',     0); 

    options.xScale = 180/pi;
    figure(4)
    printResultsToTextFiles(study,options)
end
T = readtable('../plotData/refSolutions/BCA_SHBC_MS_AS_E0_F1.txt','FileType','text', 'HeaderLines',11);
x = T.Var1;
y = T.Var2;
plot(x,y,'DisplayName', 'BEM CBM')
legend show


p = studies(1).tasks(end).task.results.abs_p;
p_ref = 10.^(y/20);
l2Error = 100*sqrt(sum((p-p_ref).^2)/sum(p_ref.^2))