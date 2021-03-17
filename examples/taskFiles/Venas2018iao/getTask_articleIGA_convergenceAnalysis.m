%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 9 and 10 in Venas2018iao
% Venas2018iao is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
misc.model = 'IL';

varCol = setIhlenburgParameters(3);
varCol{1}.meshFile = 'createNURBSmesh_EL';
parm = 1;
% misc.method = 'BA';
% misc.method = {'IENSG'};
% misc.method = {'BEM'};
% BC = {'SHBC', 'SSBC','NNBC'};
% for BC = {'SHBC', 'SSBC','NNBC'}
M_max = 5; % 7
for BC = {'SHBC'}
    misc.method = {'IE'};

%     misc.coreMethod = {'IGA'};
    misc.coreMethod = {'hp_FEM','h_FEM','C0_IGA','IGA'};
    formulation = {'BGU'};
    c_f = 1524;
    k = 1;
    omega = k*c_f;
    f = omega/(2*pi); 
    if strcmp(BC{1}, 'SHBC')
        M = 1:M_max; %1:7
    elseif strcmp(BC{1}, 'SSBC')
        M = 1:M_max-1; %1:6
    elseif strcmp(BC{1}, 'NNBC')
        M = 1:M_max-2; %1:5
    end
    alpha_s = pi;
    beta_s = 0;

    degree = 2;
    plotResultsInParaview = 0;
    calculateFarFieldPattern = 0;
    calculateVolumeError = 1;
    err.calculateSurfaceError = 1;
    prePlot.plot2Dgeometry = 0;
    prePlot.plot3Dgeometry = 0;
    loopParameters = {'M','degree','misc.method','misc.coreMethod','formulation','BC'};

    N = 6;
    
    collectIntoTasks
    
%     misc.method = {'BA'};
%     formulation = {'VL2E'};
% %     collectIntoTasks
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    misc.method = 'IE';
    misc.coreMethod = 'IGA';
    formulation = 'BGU';
    degree = 3:4;
    collectIntoTasks
    
%     misc.method = {'BA'};
%     formulation = {'VL2E'};
% %     collectIntoTasks    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    misc.method = 'IE';
    formulation = 'BGU';
    misc.coreMethod = 'linear_FEM';
    degree = 1;

    collectIntoTasks
    
%     misc.method = 'BA';
%     formulation = 'VL2E';
% %     collectIntoTasks    
end