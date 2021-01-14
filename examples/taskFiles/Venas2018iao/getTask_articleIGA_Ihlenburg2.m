function studies = getTask_articleIGA_Ihlenburg2()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Table 2 in Venas2018iao
% Venas2018iao is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

prePlot.plot2Dgeometry = 0; 
prePlot.abortAfterPlotting = 1;
scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
model = 'IL';
formulation = 'BGU';
method = 'IE';
parm = 1;
coreMethods = {'IGA','IGA','hp_FEM','linear_FEM'}; % [5, 4, 2, 1, 3]
coreMethods = {'IGA'}; % [5, 4, 2, 1, 3]
BCs = {'SHBC','SSBC','NNBC'};
BCs = {'NNBC'};
postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= 0;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'method','parm','f','formulation','M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
% postPlot(1).addCommands  = @(study,i_study,studies) addCommands_(i_study,studies);
warning('off','NURBS:weights')
M_0 = 1; % 4
for ii = 1:length(coreMethods)
    coreMethod = coreMethods(ii);
    for BC = BCs
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
        varCol{1}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, max(2^(M-1)/8-1,0)];
        if numel(varCol) > 1
            varCol{2}.refinement = @(M,t,t_fluid) [2^(M-1)-1, 2^(M-1)-1, max(round(t/t_fluid)*2^(M-1),0)];
        end
        if numel(varCol) > 2
            varCol{3}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, max(2^(M-2)-1,0)];
        end
        alpha_s = pi;
        beta_s = 0;
        plotResultsInParaview = 0;
        calculateFarFieldPattern = 1;
        calculateVolumeError = 1;
        calculateSurfaceError = 0;
        loopParameters = {'M','degree','coreMethod','BC','f'};

        collectIntoTasks
    end
end

function addCommands_(i_study,studies)
% close all
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
