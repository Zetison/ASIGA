close all
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'k',  ...
                     'yname',           'TS', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'plot', ... 
                     'lineStyle',       '-', ... 
                     'xLoopName',       'f', ...
                     'subFolderName',   '../results/Cube_P_sweep', ...
                     'legendEntries',   {{'method','M','parm','formulation'}}, ...
                     'noXLoopPrms',     1); 
    figure(2)
    printResultsToTextFiles(study,options)

    options.yname = 'error_p';
    options.axisType = 'semilogy';

    figure(4)
    printResultsToTextFiles(study,options)

    options.yname = 'surfaceError';

    figure(6)
    printResultsToTextFiles(study,options)

end

figure(6)
eigenValues = [1,2,3,4,5,6,8,9,10];
eigenValues = [eigenValues, 3,6,9];
eigenValues = pi*sqrt(eigenValues);%  Analytical eigenvalues of the interior Dirichlet/Neumann cube problem.

for i = 1:numel(eigenValues)
    yLim = ylim;
    if i < 10
        color = 'red';
    else
        color = 'blue';
    end
    plot([1,1]*eigenValues(i), yLim,'--','color',color)
end