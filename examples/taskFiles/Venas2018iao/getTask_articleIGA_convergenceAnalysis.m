function studies = getTask_articleIGA_convergenceAnalysis()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 9 and 10 in Venas2018iao
% Venas2018iao is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
misc.model = 'IL';

msh.meshFile = 'createNURBSmesh_EL';
msh.parm = 1;
% misc.method = 'BA';
% misc.method = {'IENSG'};
% misc.method = {'BEM'};
% BC = {'SHBC', 'SSBC','NNBC'};
% for BC = {'SHBC', 'SSBC','NNBC'}

postPlot(1).xname       	= 'dofs';
postPlot(1).yname        	= 'energyError';
postPlot(1).plotResults  	= 1;
postPlot(1).printResults 	= 1;
postPlot(1).axisType        = 'loglog';
postPlot(1).lineStyle   	= '*-';
postPlot(1).xLoopName     	= 'msh.M';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;

M_max = 5; % 7
for BC = {'SHBC'}
    misc.BC = BC{1};
    misc.method = {'IE'};

%     misc.coreMethod = {'IGA'};
%     misc.coreMethod = {'hp_FEM'};
    misc.coreMethod = {'hp_FEM','h_FEM','C0_IGA','IGA'};
    misc.formulation = {'BGU'};
    c_f = 1524;
    k = 1;
    misc.omega = k*c_f;
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
    ffp.alpha_s = pi;
    ffp.beta_s = 0;

    msh.degree = 2;
    para.plotResultsInParaview = 0;
    ffp.calculateFarFieldPattern = 0;
    err.calculateVolumeError = 1;
    err.calculateSurfaceError = 1;
    prePlot.plot2Dgeometry = 0;
    prePlot.plot3Dgeometry = 0;
    loopParameters = {'msh.M','msh.degree','misc.method','misc.coreMethod','misc.formulation','misc.BC'};

    iem.N = 6;
    
    collectIntoTasks
    
%     misc.method = {'BA'};
%     formulation = {'VL2E'};
% %     collectIntoTasks
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    misc.method = 'IE';
    misc.coreMethod = 'IGA';
    misc.formulation = 'BGU';
    msh.degree = 3:4;
%     collectIntoTasks
    
%     misc.method = {'BA'};
%     formulation = {'VL2E'};
% %     collectIntoTasks    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    misc.method = 'IE';
    misc.formulation = 'BGU';
    misc.coreMethod = 'linear_FEM';
    msh.degree = 1;

%     collectIntoTasks
    
%     misc.method = 'BA';
%     formulation = 'VL2E';
% %     collectIntoTasks    
end