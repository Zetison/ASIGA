function studies = getTask_Venas2018iao_Table2(M_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Table 2 in Venas2018iao
% Venas2018iao is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

if nargin < 1
    M_0 = 1; % 4
end

prePlot.plot2Dgeometry = 0; 
prePlot.plot3Dgeometry = 0; 
prePlot.abortAfterPlotting = 1;
misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
misc.model = 'IL';
misc.formulation = {'BGU'};
misc.method = {'IE'};
msh.parm = 1;
misc.coreMethods = {'IGA','IGA','hp_FEM','linear_FEM'}; % [5, 4, 2, 1, 3]
% misc.coreMethods = {'IGA'}; % [5, 4, 2, 1, 3]
% misc.coreMethods = {'hp_FEM'}; % [5, 4, 2, 1, 3]
BCs = {'SHBC','SSBC','NNBC'};
% BCs = {'NNBC'};
% BCs = {'SHBC'};
% BCs = {'SSBC'};
postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= 0;
postPlot(1).printResults 	= 0;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'msh.M';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
postPlot(1).addCommands  = @(study,i_study,studies) addCommands_(i_study,studies);

warning('off','NURBS:weights')
for i_coreM = 1:length(misc.coreMethods)
    misc.coreMethod = misc.coreMethods(i_coreM);
    for BC = BCs
        misc.BC = BC{1};
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
        misc.omega = k*c_f;

        switch misc.coreMethod{1}
            case 'IGA'
                if i_coreM == 1
                    msh.M = M_0; % 4
                    msh.degree = 3; %3
                else
                    msh.M = M_0+1;
                    msh.degree = 2:3;
                end
            case 'hp_FEM'
                msh.degree = 2;
                msh.M = M_0+1;
            case 'linear_FEM'
                msh.degree = 2;
                msh.M = M_0+2;
        end
        varCol = setIhlenburgParameters(noDomains);
        msh.meshFile = 'createNURBSmesh_EL';
        varCol{1}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, max(2^(M-1)/8-1,0)];
        if numel(varCol) > 1
            varCol{2}.refinement = @(M,t,t_fluid) [2^(M-1)-1, 2^(M-1)-1, max(round(t/t_fluid)*2^(M-1),0)];
        end
        if numel(varCol) > 2
            varCol{3}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, max(2^(M-2)-1,0)];
        end
        ffp.alpha_s = pi;
        ffp.beta_s = 0;
        
        para.plotResultsInParaview = 0;
        err.calculateFarFieldPattern = 1;
        err.calculateVolumeError = 1;
        err.calculateSurfaceError = 0;
        loopParameters = {'misc.method','misc.formulation','msh.M','msh.degree','misc.coreMethod','misc.BC','misc.omega'};

        misc.formulation = {'BGU'};
        misc.method = {'IE'};    
        collectIntoTasks
%         
%         misc.formulation = {'VL2E'};
%         misc.method = {'BA'};
%         collectIntoTasks
    end
end

function addCommands_(i_study,studies)
printForLaTeX = 0;
for j = 1:3
    tasks = studies(9+j).tasks;
    tasks(end+1) = studies(6+j).tasks;
    tasks(end+1) = studies(3+j).tasks(1);
    tasks(end+1) = studies(j).tasks(1);
    tasks(end+1) = studies(3+j).tasks(2);
    studies2(j).tasks = tasks;
end
if ~printForLaTeX
    fprintf('%-20s\t  %6s\t  %6s \t  %5s  %6s  %5s\n', ...
             ' ', 'noElems', 'dofs', 't_sys', 't_sol', 'energyError');
end
for j = 1:3
    tasks = studies2(j).tasks;
    for ii = 1:numel(tasks)
        task = tasks(ii).task;
        M = task.msh.M;
        degree = task.msh.degree;
        t_sys = task.varCol{1}.timeBuildSystem;
        t_sol = task.varCol{1}.tot_time - t_sys;
        noElems = task.varCol{1}.totNoElems;
        dofs = task.dofs;
        switch task.misc.coreMethod
            case 'IGA'
                if printForLaTeX
                    fprintf('Mesh ${\\cal M}_{%d,%d,%d}^{\\textsc{iga}}$\t\t & %6d\t & %6d \t & %5.2f & %6.2f & %5.2f \\cr\n', ...
                             M, degree, degree-1, noElems, dofs, t_sys, t_sol, task.results.energyError);
                else
                    fprintf('Mesh M_{%d,%d,%d}^iga\t  %6d\t  %6d \t  %5.2f  %6.2f  %5.2f\n', ...
                             M, degree, degree-1, noElems, dofs, t_sys, t_sol, task.results.energyError);
                end
            otherwise
                if printForLaTeX
                    fprintf('Mesh ${\\cal M}_{%d,%d,\\mathrm{i}}^{\\textsc{fem}}$\t & %6d\t & %6d \t & %5.2f & %6.2f & %5.2f \\cr\n', ...
                             M, degree, noElems, dofs, t_sys, t_sol, task.results.energyError);
                else
                    fprintf('Mesh M_{%d,%d,i}^fem\t  %6d\t  %6d \t  %5.2f  %6.2f  %5.2f\n', ...
                             M, degree, noElems, dofs, t_sys, t_sol, task.results.energyError);
                end
        end
    end
    fprintf('\n')
end
