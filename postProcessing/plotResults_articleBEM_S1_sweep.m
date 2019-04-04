
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'k',  ...
                     'yname',           'TS', ...
                     'plotResults', 	1, ... 
                     'printResults',	1, ... 
                     'axisType',        'plot', ... 
                     'lineStyle',       '-', ... 
                     'xLoopName',       'f', ...
                     'subFolderName',   '../results/S1_sweep_BEM', ...
                     'legendEntries',   {{'method','M','parm','formulation'}}, ...
                     'noXLoopPrms',     any(strcmp(study.loopParameters,'f'))); 
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
eigenValues = [  pi  % Analytical eigenvalues of the interior Dirichlet sphere example.
               2*pi
               3*pi
               4.493409457909065
               7.725251836937707
               5.763459196894550
               9.095011330476353
               6.987932000500519
               8.182561452571243
               9.355812111042747
               4.493409457909064 % Analytical eigenvalues of the interior Neumann sphere example.
               7.725251836937708
               2.081575977818101
               5.940369990572713
               9.205840142936665
               3.342093657365695
               7.289932304093351
               4.514099647032276
               8.583754956365768
               5.646703620436797
               9.840446043040137               
               6.756456330204129
               7.851077679474405
               8.934838878352839
               ]';
for i = 1:numel(eigenValues)
    yLim = ylim;
    if i < 10
        color = 'red';
    else
        color = 'blue';
    end
    plot([1,1]*eigenValues(i), yLim,'--','color',color)
end