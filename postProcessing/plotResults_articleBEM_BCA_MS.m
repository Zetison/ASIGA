
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'alpha',  ...
                     'yname',           'TS', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'plot', ... 
                     'lineStyle',       '-', ... 
                     'subFolderName',   '../results/BCA_MS', ...
                     'legendEntries',   {{'method','coreMethod','formulation','M','degree','f'}}, ...
                     'noXLoopPrms',     0); 

    options.xScale = 180/pi;
    figure(4)
    printResultsToTextFiles(study,options)
end
for res = 1:5
    T = readtable(['../plotData/refSolutions/WTD71/Mo_RES' num2str(res) '_1000Hz.txt'],'FileType','text', 'HeaderLines',0);
    x = T.Var1;
    y = T.Var2;
%     x(idx) = [];
%     y(idx) = [];
    figure(4)
    plot(x,y,'DisplayName',['WTD71, M = ' num2str(res) ', f = ' num2str(f)])
    legend('off');
    legend('show');
    hold on
end
T = readtable('../plotData/refSolutions/WTD71/BC_HWBC_MS_AS_E0_F1.txt','FileType','text', 'HeaderLines',7);
x = T.Var1;
y = T.Var2;
plot(x,y,'DisplayName','WTD old')
T = readtable('../plotData/refSolutions/BC_SHBC_MS_AS_E0_F1.txt','FileType','text', 'HeaderLines',7);
x = T.Var1;
y = T.Var2;
plot(x,y,'DisplayName','WTD old 2')