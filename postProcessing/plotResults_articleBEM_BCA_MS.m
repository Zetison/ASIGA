close all
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'alpha',  ...
                     'yname',           'TS', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'plot', ... 
                     'lineStyle',       '-', ... 
                     'subFolderName',   '../results/articleBEM_BCA_MS', ...
                     'legendEntries',   {{'method','coreMethod','formulation','M','degree','f'}}, ...
                     'noXLoopPrms',     0); 

    options.xScale = 180/pi;
    figure(4)
    printResultsToTextFiles(study,options)
end
f = 1000;
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

p_ref6 = studies(1).tasks(1).task.results.abs_p;
p_ref5 = studies(2).tasks(end).task.results.abs_p;
T = readtable(['../plotData/refSolutions/WTD71/Mo_RES' num2str(5) '_1000Hz.txt'],'FileType','text', 'HeaderLines',0);

y = T.Var2;
x = T.Var1;
l2Error = 100*sqrt(sum((p_ref5-p_ref6).^2)/sum(p_ref6.^2))
l2ErrorWTD = 100*sqrt(sum((10.^(y(3601:end)/20)-p_ref6).^2)./sum(p_ref6.^2))

