close all
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'k',  ...
                     'yname',           'TS', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'semilogx', ... 
                     'lineStyle',       '-', ... 
                     'xLoopName',       'f', ...
                     'subFolderName',   '../results/Fillinger', ...
                     'legendEntries',   {{'method','M','parm','coreMethod'}}, ...
                     'noXLoopPrms',     0); 

    figure(2)
    printResultsToTextFiles(study,options)

    options.yname = 'error_p';
    options.axisType = 'loglog';

    figure(4)
    printResultsToTextFiles(study,options)

    options.yname = 'error_pAbs';
    options.axisType = 'loglog';

    figure(5)
    printResultsToTextFiles(study,options)
end
nFreqs = 1000;
k = 10.^linspace(-1,4,nFreqs).';
P_inc = 1;
R = 1;
p_ref = exactKDT(k,P_inc,R);
printResultsToFile2('../results/Fillinger/kirchhoff', k, 20*log10(p_ref), 'k', 'TS');
p_ref2 = P_inc*R/2*exp(-2*1i*k*R);
printResultsToFile2('../results/Fillinger/asymptotic', k, 20*log10(p_ref2), 'k', 'TS');
% p3 = readLaTeXFormat('results/_studies/Fillinger/S1_KDT_M6_linear_FEM_error_pVSk.txt');
for i = 1:numel(studies.tasks)
    p = studies.tasks(i).task.results.p.';
    M = studies.tasks(i).task.M;
    E = 100*abs(p_ref-p)./abs(p_ref);
    figure(45)
    loglog(k,E)
    hold on
    printResultsToFile2(['../results/Fillinger/S1_KDT_M' num2str(M) '_parm2_linear_FEM_KirchRef'], k, E, 'k', 'error_p');
end


c_f = 1500;
omega = k*c_f;
options = struct('d_vec', -[0,0,1]',... 
                 'omega', omega, ...
                 'P_inc', P_inc, ...
                 'rho_f', 1, ...
                 'calc_farField', 1, ...
                 'c_f', c_f);
             
             
v = R(1)*[0,0,1];
data = e3Dss(v, options);
figure(2)

E = 100*abs(data.p.'-p_ref)./abs(data.p).';
E2 = 100*abs(data.p.'-p_ref2)./abs(data.p).';
E2abs = 100*abs(abs(data.p.')-abs(p_ref2))./abs(data.p).';
figure(46)
loglog(k,E,k,E2)
printResultsToFile2('../results/Fillinger/kirchhoff', k, E, 'k', 'error_p');
printResultsToFile2('../results/Fillinger/asymptotic', k, E2, 'k', 'error_p');
printResultsToFile2('../results/Fillinger/asymptotic', k, E2abs, 'k', 'error_pAbs');

