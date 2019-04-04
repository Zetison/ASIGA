
figure(2)
xlim([6, 54])
ylim([5e-3,2e2])
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'dofsAlg',  ...
                     'yname',           'H1sError', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'semilogy', ... 
                     'lineStyle',       '*-', ... 
                     'xLoopName',       'degree', ...
                     'subFolderName',   '../results/S1_Demkowicz', ...
                     'noXLoopPrms',     1); 

    figure(2)
    printResultsToTextFiles(study,options)
% 
%     figure(3)
%     options.yname = 'L2Error';
%     printResultsToTextFiles(study,options)

%     figure(4)
%     options.xname = 'alpha';
%     options.yname = 'error_p';
%     options.noXLoopPrms = 0;
%     options.lineStyle = '-';
%     printResultsToTextFiles(study,options)
end
figure(2)
xlim([6, 54])
ylim([5e-3,2e2])

Michler_Inf = importdata('e3Dss/models/Michler2007itp/Figure6b_inf.csv');
Michler_PML = importdata('e3Dss/models/Michler2007itp/Figure6b_PML.csv');

plot(Michler_Inf(:,1),Michler_Inf(:,2),'*-','DisplayName','Demkowicz infinite elements')
plot(Michler_PML(:,1),Michler_PML(:,2),'*-','DisplayName','Demkowicz PML')


xname = options.xname;
yname = options.yname;
model = study.tasks(1).task.model;
savefig([options.subFolderName '/_' model '_' yname 'VS' xname])