
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'alpha',  ...
                     'yname',           'TS', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'plot', ... 
                     'lineStyle',       '-', ... 
                     'subFolderName',   '../results/BCA', ...
                     'legendEntries',   {{'method','coreMethod','formulation','M','degree','f'}}, ...
                     'noXLoopPrms',     0); 

    options.xScale = 180/pi;
    figure(4)
    printResultsToTextFiles(study,options)
end
alpha = studies(1).tasks.task.alpha;
alpha = alpha(1:end-1);
% for res = [4,8]
l2errorCOMSOL = zeros(1,4);
p_ref = studies(1).tasks.task.results.abs_p;
counter = 1;
for res = [10,20,40,80]
%     f = 1000;
    f = 100;
    alpha_s = 240;
    T = readtable(['../../comsol/models/BC/BeTSSi_mod/BC_resolution_' num2str(res) '_noPMLayers_' num2str(10) '_A' num2str(alpha_s) '_F' num2str(f) '.txt'],'FileType','text', 'HeaderLines',8);
%             T = readtable(['../comsol/models/BC/BeTSSi_mod/BETSSI~' num2str(res) '.TXT'],'FileType','text', 'HeaderLines',8);
    x = T.Var1;
    y = T.Var2;
    figure(1)
    plot(x,y,'DisplayName',['COMSOL, M = ' num2str(res) ', f = ' num2str(f)])
    legend('off');
    legend('show');
    hold on
%     Error = 100*abs(10.^(y/20)-p_ref(1:end-1))./p_ref(1:end-1);
    Error = 100*abs(10.^(y/20)-p_ref(1:end-1))./max(p_ref(1:end-1));
    l2errorCOMSOL(counter) = 100*sqrt(sum((10.^(y/20)-p_ref(1:end-1)).^2)./sum(p_ref(1:end-1).^2));
    counter = counter + 1;
    printResultsToFile2(['../results/BCA/COMSOL_Error_res' num2str(res)], 180/pi*alpha.', Error)
    figure(42)
    semilogy(180/pi*alpha,Error,'DisplayName',['COMSOL res' num2str(res)])
    hold on
    
end
l2errorCOMSOL %  1.6563    0.3034    0.2532    0.2420
savefig([options.subFolderName '/_comparison'])
figure(42)
p_ref = studies(1).tasks.task.results.abs_p;
alpha = studies(1).tasks.task.alpha;
l2errorIGA = zeros(1,6);
for task_i = 1:6 
    p = studies(2).tasks(task_i).task.results.abs_p;
    degree = studies(2).tasks(task_i).task.degree;
    M = studies(2).tasks(task_i).task.M;
    Error = 100*abs(p-p_ref)./max(p_ref);
    l2errorIGA(task_i) = 100*sqrt(sum((p-p_ref).^2)./sum(p_ref.^2));
    filename = ['M' num2str(M) 'degree' num2str(degree) 'f100'];
    printResultsToFile2(['../results/BCA/BCA_BEM_IGA_CCBIE_' filename], 180/pi*alpha.', Error)
    semilogy(180/pi*alpha,Error,'DisplayName',filename)
    hold on
end
l2errorIGA % 0.1226    0.0412    0.0600    0.0206    0.0347    0.0089
yLabel = '$$\frac{||p|-|p_h||}{|p|} [\%]$$';
ylabel(yLabel,'interpreter','latex')
xlabel('$$\alpha$$','interpreter','latex')
legend show
xlim([0,360])
savefig([options.subFolderName '/_error'])

fprintf('\n\n')
for i = 1:6
    fprintf('\t\t${\\cal M}_{%d,%d,%d}^{\\textsc{igabem}}$ & %d & %d & %g & %g & %g & %d\\\\\n', studies(2).tasks(i).task.M, studies(2).tasks(i).task.degree, studies(2).tasks(i).task.degree-1, studies(2).tasks(i).task.varCol.noElems, studies(2).tasks(i).task.varCol.dofs, studies(2).tasks(i).task.varCol.h_max, studies(2).tasks(i).task.varCol.nepw, l2errorIGA(i), round(studies(2).tasks(i).task.varCol.tot_time)); 
end
i = 1;
fprintf('\t\t${\\cal M}_{%d,%d,%d}^{\\textsc{igabem}}$ & %d & %d & %g & %g & 0 & %d\\\\\n', studies(1).tasks(i).task.M, studies(1).tasks(i).task.degree, studies(1).tasks(i).task.degree-1, studies(1).tasks(i).task.varCol.noElems, studies(1).tasks(i).task.varCol.dofs, studies(1).tasks(i).task.varCol.h_max, studies(1).tasks(i).task.varCol.nepw, round(studies(1).tasks(i).task.varCol.tot_time)); 
