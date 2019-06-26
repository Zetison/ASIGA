

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'IL';

parm = 1;
% method = 'BA';
% method = {'IENSG'};
% method = {'BEM'};
% BC = {'SHBC', 'SSBC','NNBC'};
for BC = {'SHBC'} %{'SHBC', 'SSBC','NNBC'}
    method = {'IE'};

%     coreMethod = {'C0_IGA','hp_FEM'};
%     coreMethod = {'IGA'};
    coreMethod = {'hp_FEM','h_FEM','C0_IGA','IGA'};
    formulation = {'BGU'};
%     formulation = 'VL2E';
%     formulation = 'SL2E';
    c_f = 1524;
    k = 1;
    omega = k*c_f;
    f = omega/(2*pi); 
    if strcmp(BC{1}, 'SHBC')
        M = 1:7; %3:6
%         M = 1:2; %3:6
    elseif strcmp(BC{1}, 'SSBC')
        M = 1:6; %3:6
    elseif strcmp(BC{1}, 'NNBC')
        M = 1:5; %3:6
    end
%     M = 7;
    alpha_s = pi;
    beta_s = 0;

    degree = 2;
    plotResultsInParaview = 0;
    calculateFarFieldPattern = 0;
    calculateVolumeError = 1;
    calculateSurfaceError = 1;
    plot2Dgeometry = 0;
    plot3Dgeometry = 0;
    loopParameters = {'M','degree','method','coreMethod','formulation','BC'};

    N = 5;
    
    collectIntoTasks
    
    method = {'BA'};
    formulation = 'VL2E';
    collectIntoTasks
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    method = {'IE'};
    coreMethod = {'IGA'};
    formulation = {'BGU'};
    degree = 3:4;
%     degree = 4;
    collectIntoTasks
    
    method = {'BA'};
    formulation = {'VL2E'};
    collectIntoTasks    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    method = {'IE'};
    formulation = {'BGU'};
    coreMethod = {'linear_FEM'};
    degree = 1;

    collectIntoTasks
    
    method = {'BA'};
    formulation = {'VL2E'};
    collectIntoTasks    
end