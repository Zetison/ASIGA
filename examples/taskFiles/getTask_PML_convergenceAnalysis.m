function studies = getTask_PML_convergenceAnalysis()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 9 and 10 in Venas2018iao
% Venas2018iao is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
misc.model = 'IL';
noCoresToUse = 12;

msh.meshFile = 'createNURBSmesh_EL';
msh.parm = 1;
% misc.method = 'BA';
% misc.method = {'IENSG'};
% misc.method = {'BEM'};
% BC = {'SHBC', 'SSBC','NNBC'};
% for BC = {'SHBC', 'SSBC','NNBC'}

prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;
prePlot.plotControlPolygon  = 1;       % Plot the control polygon for the NURBS mesh
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
% prePlot.colorFun = @(v) abs(norm2(v)-(r_a+t_PML));
prePlot.resolution = [20,20,0];
warning('off','NURBS:weights')

% postPlot(1).xname       	= 'nepw';
postPlot(1).xname       	= 'dofs';
postPlot(1).yname        	= 'energyError';
postPlot(1).plotResults  	= 1;
postPlot(1).printResults 	= 1;
postPlot(1).axisType        = 'loglog';
postPlot(1).lineStyle   	= '*-';
postPlot(1).xLoopName     	= 'msh.M';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;


msh.meshFile = 'createNURBSmesh_EL';
msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
msh.refineThetaOnly = false;

for method = {'PML'} %{'IE','PML'}
    misc.method = method{1};

    if strcmp(method{1},'PML')
        misc.formulation = {'GSB'};
    else
        misc.formulation = {'BGU'};
    end

    M_max = 6; % 7
    for BC = {'SHBC'}
        misc.BC = BC{1};
        misc.coreMethod = {'IGA'};
%         misc.coreMethod = {'hp_FEM'};
        misc.coreMethod = {'hp_FEM','h_FEM','C0_IGA','IGA'};
%         misc.coreMethod = {'C0_IGA','IGA'};
        c_f = 1524;
        k = 1;
        misc.omega = k*c_f;
        msh.degree = 2;
        if strcmp(BC{1}, 'SHBC')
            varCol = setIhlenburgParameters(1);
            msh.M = 1:M_max; %1:7
        elseif strcmp(BC{1}, 'SSBC')
            varCol = setIhlenburgParameters(2);
            msh.M = 1:M_max-1; %1:6
        elseif strcmp(BC{1}, 'NNBC')
            varCol = setIhlenburgParameters(3);
            msh.M = 1:M_max-2; %1:5
        end
%         msh.M = (M_max-1):M_max; %1:5
%         msh.M = 5;
        iem.N = 6;
        pml.t = 0.25*varCol{1}.R_i;         % thickness of PML
        pml.dirichlet = true;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
        pml.sigmaType = 3;   	% sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n, sigmaType = 3: sigma(xi) = C/(1-xi)^n
        if pml.sigmaType == 3
            pml.n = 1;            	% polynomial order
        else
            pml.n = 2;            	% polynomial order
        end
        misc.r_a = 1.25*varCol{1}.R_i;
%         varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-4)-1, max(iem.N - msh.degree,2^(M-3)-1)];
%         varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-4)-1, 2^(M-4)-1];
        varCol{1}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, 2^(M-4)-1, 2^(M-4)-1];
        ffp.alpha_s = 0;
        ffp.beta_s = -pi/2;

        para.plotResultsInParaview = 0;
        ffp.calculateFarFieldPattern = 0;
        err.calculateVolumeError = 1;
        err.calculateSurfaceError = 1;
        loopParameters = {'msh.M','msh.degree','misc.method','misc.coreMethod','misc.formulation','misc.BC'};


        collectIntoTasks

    %     misc.method = {'BA'};
    %     formulation = {'VL2E'};
    % %     collectIntoTasks

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        misc.coreMethod = 'IGA';
        msh.degree = 3:4;
    %     collectIntoTasks

    %     misc.method = {'BA'};
    %     formulation = {'VL2E'};
    % %     collectIntoTasks    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        misc.coreMethod = 'linear_FEM';
        msh.degree = 1;

        collectIntoTasks

    %     misc.method = 'BA';
    %     formulation = 'VL2E';
    % %     collectIntoTasks    
    end
end