close all
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'alpha',  ...
                     'yname',           'TS', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'plot', ... 
                     'lineStyle',       '-', ... 
                     'xLoopName',       'N', ...
                     'subFolderName',   '../results/PhD_M3/', ...
                     'legendEntries',   {{'method','coreMethod','formulation','M','N'}}, ...
                     'noXLoopPrms',     0); 

    options.xScale = 180/pi;
    figure(5)
    printResultsToTextFiles(study,options)
end
% if 1
%     f = 1000;
%     switch f
%         case 100
%             files = dir('../plotData/M3_HWBC_MS_old/f_100_*.dat');
%             for file = files'
% %                 temp = importdata(['../plotData/M3_HWBC_MS_old/' file.name], ' ', 1);
%                 temp = readtable(['../plotData/M3_HWBC_MS_old/' file.name],'FileType','text', 'HeaderLines',1);
%                 x = temp.Var1;
%                 y = temp.Var2;
%                 plot(x,y,'DisplayName',file.name(1:end-4))
%                 legend('off');
%                 l = legend('show');
%                 hold on
%             end
%             set(l, 'Interpreter', 'none')
%         case 1000
%             files = dir('../plotData/refSolutions/M3_HWBC_MS_0_1*.txt');
%             for file = files'
%                 temp = importdata(['../plotData/refSolutions/' file.name], ',', 7);
%             %             T = readtable('plotData/refSolutions/PH_HWBC_MS_AS_E0_F1.txt','FileType','text', 'HeaderLines',12);
%                 x = temp.data(:,1);
%                 y = temp.data(:,2);
%                 plot(x,y,'DisplayName',file.name(1:end-4))
%                 legend('off');
%                 l = legend('show');
%                 hold on
%             end
%             set(l, 'Interpreter', 'none')
%     end
% else
%     
%     T = readtable('plotData/M3_HWBC_MS_old/f_100_fullRange_mesh05.dat','FileType','text', 'HeaderLines',1);
%     T2 = readtable('plotData/M3_HWBC_MS_old/f_100_fullRange_mesh2.dat','FileType','text', 'HeaderLines',1);
%     hold on 
%     plot(T.Var1,T.Var2,'DisplayName','IGA IE M=5')
%     plot(T2.Var1,T2.Var2,'DisplayName','IGA IE M=2')
%     legend('off');
%     legend('show');
% end
% 
% figure(42)
% p_ref = studies(1).tasks(6).task.results.abs_p;
% alpha = studies(1).tasks(6).task.alpha;
% l2errorIGA = zeros(1,6);
% l2errorIGA2 = zeros(1,6);
% for task_i = 1:6 
%     p = studies(1).tasks(task_i).task.results.abs_p;
%     degree = studies(1).tasks(task_i).task.degree;
%     formulation = studies(1).tasks(task_i).task.formulation;
%     M = studies(1).tasks(task_i).task.M;
%     Error = 100*abs(p-p_ref)./max(p_ref);
%     l2errorIGA(task_i) = 100*sqrt(sum((p-p_ref).^2)./sum(p_ref.^2));
% %     l2errorIGA2(task_i) = 100*sqrt(sum((p(1:end-1)-p_ref2).^2)./sum(p_ref2.^2));
%     filename = ['M' num2str(M) 'degree' num2str(degree) 'f' num2str(f) formulation];
% %     printResultsToFile2(['../results/articleBEM_BCA/BCA_BEM_IGA_CCBIE_' filename], 180/pi*alpha.', Error)
%     semilogy(180/pi*alpha,Error,'DisplayName',filename)
%     hold on
% end
% l2errorIGA 
% yLabel = '$$\frac{||p|-|p_h||}{|p|} [\%]$$';
% ylabel(yLabel,'interpreter','latex')
% xlabel('$$\alpha$$','interpreter','latex')
% legend show
% xlim([0,180])
% savefig([options.subFolderName '/_error'])


