
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'alpha',  ...
                     'yname',           'TS', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'plot', ... 
                     'lineStyle',       '-', ... 
                     'xLoopName',       'N', ...
                     'subFolderName',   'results/M3/', ...
                     'legendEntries',   {{'method','coreMethod','formulation','M','degreeElev','f','N'}}, ...
                     'noXLoopPrms',     0); 

    options.xScale = 180/pi;
    figure(5)
    printResultsToTextFiles(study,options)
end
if 0
    files = dir('plotData/refSolutions/M3_HWBC_MS_0_1*.txt');
    for file = files'
        temp = importdata(['plotData/refSolutions/' file.name], ',', 7);
    %             T = readtable('plotData/refSolutions/PH_HWBC_MS_AS_E0_F1.txt','FileType','text', 'HeaderLines',12);
        x = temp.data(:,1);
        y = temp.data(:,2);
        plot(x,y,'DisplayName',file.name(1:end-4))
        legend('off');
        l = legend('show');
        hold on
    end
    set(l, 'Interpreter', 'none')
else
    
    T = readtable('plotData/M3_HWBC_MS_old/f_100_fullRange_mesh05.dat','FileType','text', 'HeaderLines',1);
    T2 = readtable('plotData/M3_HWBC_MS_old/f_100_fullRange_mesh2.dat','FileType','text', 'HeaderLines',1);
    hold on 
    plot(T.Var1,T.Var2,'DisplayName','IGA IE M=5')
    plot(T2.Var1,T2.Var2,'DisplayName','IGA IE M=2')
    legend('off');
    legend('show');
end