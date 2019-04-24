
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'alpha',  ...
                     'yname',           'TS', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'plot', ... 
                     'lineStyle',       '-', ... 
                     'xLoopName',       'N', ...
                     'subFolderName',   'results/BC/allStudies', ...
                     'legendEntries',   {{'method','coreMethod','formulation','M','degreeElev','f','N'}}, ...
                     'noXLoopPrms',     0); 

    options.xScale = 180/pi;
    figure(5)
    printResultsToTextFiles(study,options)
end
for largerDomain = 0 %[0,1]
    for res = [4,8,16,32]
        if res > 10 && largerDomain
            continue
        end
        f = 100;
        if largerDomain
%                     T = readtable(['../comsol/models/BC/BC_F' num2str(f) '_resolution_' num2str(res) '_largeDomain.txt'],'FileType','text', 'HeaderLines',1);
            T = readtable(['../comsol/models/BC/BC_F' num2str(f) '_resolution_' num2str(res) '_old.txt'],'FileType','text', 'HeaderLines',8);
        else
            T = readtable(['../comsol/models/BC/BC_F' num2str(f) '_resolution_' num2str(res) '.txt'],'FileType','text', 'HeaderLines',8);
        end
        x = T.Var1;
        y = T.Var2;
        plot(x,y,'DisplayName',['COMSOL, M = ' num2str(res) ', f = 100, largerDomain = ' num2str(largerDomain)])
        legend('off');
        legend('show');
        hold on
    end
end
for res = [5,10]
    f = 500;
    T = readtable(['../comsol/models/BC/BC_F' num2str(f) '_resolution_' num2str(res) '.txt'],'FileType','text', 'HeaderLines',8);
    x = T.Var1;
    y = T.Var2;
    plot(x,y,'DisplayName',['COMSOL, M = ' num2str(res) ', f = 500, largerDomain = ' num2str(largerDomain)])
    legend('off');
    legend('show');
    hold on
end
%         for res = [10,20]
%             f = 100;
%             T = readtable(['../comsol/models/BC/BC_F' num2str(f) '_resolution_' num2str(res) '_old.txt'],'FileType','text', 'HeaderLines',8);
%             x = T.Var1;
%             y = T.Var2;
%             plot(x,y,'DisplayName',['COMSOL, M = ' num2str(res) ', f = 100, largerDomain = 1'])
%             legend('off');
%             legend('show');
%             hold on
%         end
for res = [10,20]
    f = 100;
    T = readtable(['../comsol/models/BC/BC_F' num2str(f) '_resolution_' num2str(res) '.TXT'],'FileType','text', 'HeaderLines',8);
    x = T.Var1;
    y = T.Var2;
    plot(x,y,'DisplayName',['COMSOL, M = ' num2str(res) ', f = 100, free tetrahedral mesh'])
    legend('off');
    legend('show');
    hold on
end

%         figure(6)
% 
%         options.noXLoopPrms = 0;
%         options.lineStyle = '-';
%         options.xname = 'alpha';
%         options.yname = 'error_pAbs';
%         options.axisType = 'semilogy';
%         options.xScale = 180/pi;
%         for study_i = 1:numel(studies)  
%             study = studies(study_i);
% 
%             options.xScale = 180/pi;
%             printResultsToTextFiles(study,options)
%         end

%         figure(7)
% 
%         options.noXLoopPrms = 0;
%         options.lineStyle = '-';
%         options.xname = 'alpha';
%         options.yname = 'error_p';
%         options.axisType = 'semilogy';
%         options.xScale = 180/pi;
%         for study_i = 1:numel(studies)  
%             study = studies(study_i);
% 
%             options.xScale = 180/pi;
%             printResultsToTextFiles(study,options)
%         end