%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Table 2 in Venas2018iao
% Venas2018iao is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
model = 'IL';
formulation = 'BGU';
method = 'IE';
parm = 1;
coreMethods = {'IGA','IGA','hp_FEM','linear_FEM'}; % [5, 4, 2, 1, 3]
% coreMethods = {'linear_FEM'}; % [5, 4, 2, 1, 3]
postPlot(1).plotResults = false;
postPlot(1).printResults = false;
postPlot(1).addCommands  = @(study,i_study,studies) addCommands_(i_study,studies);
M_0 = 1; % 4
for ii = 1:length(coreMethods) %{'IGA'}
    coreMethod = coreMethods(ii);
    for BC = {'SHBC','SSBC','NNBC'} %
        switch BC{1}
            case 'SHBC'
                noDomains = 1;
            case 'SSBC'
                noDomains = 2;
            case 'NNBC'
                noDomains = 3;
        end
        k = 1;
        c_f = 1524;
        omega = k*c_f;
        f = omega/(2*pi); 

        switch coreMethod{1}
            case 'IGA'
                if ii == 1
                    M = M_0; % 4
                    degree = 3; %3
                else
                    M = M_0+1;
                    degree = 2:3;
                end
            case 'hp_FEM'
                degree = 2;
                M = M_0+1;
            case 'linear_FEM'
                degree = 2;
                M = M_0+2;
        end
        varCol = setIhlenburgParameters(noDomains);
        varCol{1}.meshFile = 'createNURBSmesh_EL';
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

function addCommands_(i_study,studies)
if i_study == 1
    for j = 1:3
        tasks = studies(9+j).tasks;
        tasks(end+1) = studies(6+j).tasks;
        tasks(end+1) = studies(3+j).tasks(1);
        tasks(end+1) = studies(j).tasks(1);
        tasks(end+1) = studies(3+j).tasks(2);
        studies2(j).tasks = tasks;
    end
    for j = 1:3
        tasks = studies2(j).tasks;
        for ii = 1:numel(tasks)
            task = tasks(ii).task;
            M = task.M;
            degree = task.degree;
            t_sys = task.varCol{1}.timeBuildSystem;
            t_sol = task.varCol{1}.tot_time - t_sys;
            noElems = task.varCol{1}.totNoElems;
            dofs = task.varCol{1}.dofs;
            switch task.coreMethod
                case 'IGA'
                    fprintf('Mesh ${\\cal M}_{%d,%d,%d}^{\\textsc{iga}}$\t\t\t & %6d\t & %6d \t & %5.2f & %6.2f & %5.2f \\cr\n', ...
                             M, degree, degree-1, noElems, dofs, t_sys, t_sol, task.results.energyError);
                otherwise
                    fprintf('Mesh ${\\cal M}_{%d,%d,\\mathrm{i}}^{\\textsc{fem}}$\t & %6d\t & %6d \t & %5.2f & %6.2f & %5.2f \\cr\n', ...
                             M, degree, noElems, dofs, t_sys, t_sol, task.results.energyError);
            end
        end
        fprintf('\n')
    end
end
end
