scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
model = 'IL';
formulation = 'BGU';
method = 'IE';
parm = 1;
coreMethods = {'IGA','IGA','hp_FEM','linear_FEM'}; % [5, 4, 2, 1, 3]
% coreMethods = {'IGA'}; % [5, 4, 2, 1, 3]
for i_coreM = 1:length(coreMethods) %{'IGA'}
    coreMethod = {coreMethods{i_coreM}};
    for BC = {'SHBC','SSBC','NNBC'} %
        k = 1;
        c_f = 1524;
        omega = k*c_f;
        f = omega/(2*pi); 

        switch coreMethod{1}
            case 'IGA'
                if i_coreM == 1
                    M = 4; % 4
                    degree = 3; %3
                else
                    M = 5;
                    degree = 2:3;
                end
            case 'hp_FEM'
                degree = 2;
                M = 5;
            case 'linear_FEM'
                degree = 2;
                M = 6;
        end
        alpha_s = pi;
        beta_s = 0;
        alpha = alpha_s;
        beta = beta_s;
        plotResultsInParaview = 0;
        calculateFarFieldPattern = 1;
        calculateVolumeError = 1;
        calculateSurfaceError = 0;
        loopParameters = {'M','degree','coreMethod','BC','f'};

        collectIntoTasks
    end
end
return

tasks = studies(1).tasks;
tasks(end+1:end+2) = studies(2).tasks;
tasks(end+1) = studies(3).tasks;
tasks(end+1) = studies(4).tasks;
for i = 1:numel(tasks)
    varCol = tasks(i).task.varCol;
    M = varCol.M;
    degree = varCol.degree;
    t_sys = varCol.timeBuildSystem;
    t_sol = varCol.tot_time - t_sys;
    switch varCol.coreMethod
        case 'IGA'
            fprintf('Mesh ${\\cal M}_{%d,%d,\\mathrm{i}}^{\\textsc{fem}}$\t & %d\t & %d \t & %5.2f & %5.2f & %5.2f \\cr', ...
                     M, degree, degree-1,varCol.noElems, varCol.actualNoDofs, varCol.noElems, t_sys, t_sol, varCol.results.energyError);
        otherwise
            fprintf('Mesh ${\\cal M}_{%d,%d,%d}^{\\textsc{iga}}$\t & %d\t & %d \t & %5.2f & %5.2f & %5.2f \\cr', ...
                     M, degree, varCol.noElems, varCol.actualNoDofs, varCol.noElems, t_sys, t_sol, varCol.results.energyError);
    end
end
