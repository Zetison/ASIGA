

saveStudies = false;

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
misc.model = 'IL';

parm = 1;
prePlot.plot2Dgeometry = 0;
prePlot.plot3Dgeometry = 0;
if 1
    misc.method = {'BEM'};
    formulation = 'GCBIE';
    solveForPtot = true;
else
    misc.method = {'IE'};
    formulation = 'BGU';
    solveForPtot = false;
    calculateVolumeError = 1;
    N = 6;
end
clearGlobalMatrices = false;
BC = 'NNBC';
% BC = 'SSBC';
c_f = 1524;
k = 0.1; % k = 1
omega = k*c_f;
f = omega/(2*pi); 
M = 1:3; %1:5
% M = 1; %1:5
alpha_s = pi;
beta_s = 0;
err.calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
plotResultsInParaview = 0;

extraGP = 0; % extra quadrature points
degree = 4;
loopParameters = {'M','misc.method'};
collectIntoTasks
    
misc.method = {'BA'};
% formulation = 'SL2E';
formulation = 'VL2E';
collectIntoTasks

% misc.method = 'BA';
% misc.method = {'IENSG'};
% misc.method = {'BEM'};
% BC = {'SHBC', 'SSBC','NNBC'};
% for BC = {'SHBC', 'SSBC','NNBC'}
%     misc.method = {'IE'};
% 
% %     misc.coreMethod = {'C0_IGA','hp_FEM'};
% %     misc.coreMethod = {'IGA'};
%     misc.coreMethod = {'hp_FEM','h_FEM','C0_IGA','IGA'};
%     formulation = {'BGU'};
% %     formulation = 'VL2E';
% %     formulation = 'SL2E';
%     c_f = 1524;
%     k = 1;
%     omega = k*c_f;
%     f = omega/(2*pi); 
%     if strcmp(BC{1}, 'SHBC')
%         M = 1:7; %1:7
%     elseif strcmp(BC{1}, 'SSBC')
%         M = 1:6; %1:6
%     elseif strcmp(BC{1}, 'NNBC')
%         M = 1:5; %1:5
%     end
% %     M = 1:3;
%     alpha_s = pi;
%     beta_s = 0;
% 
%     degree = 2;
%     plotResultsInParaview = 0;
%     calculateFarFieldPattern = 0;
%     calculateVolumeError = 1;
%     err.calculateSurfaceError = 1;
%     prePlot.plot2Dgeometry = 0;
%     prePlot.plot3Dgeometry = 0;
%     loopParameters = {'M','degree','misc.method','misc.coreMethod','formulation','BC'};
% 
%     N = 6;
%     
%     collectIntoTasks
%     
%     misc.method = {'BA'};
%     formulation = {'VL2E'};
%     collectIntoTasks
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     misc.method = {'IE'};
%     misc.coreMethod = {'IGA'};
%     formulation = {'BGU'};
%     degree = 3:4;
% %     degree = 4;
%     collectIntoTasks
%     
%     misc.method = {'BA'};
%     formulation = {'VL2E'};
%     collectIntoTasks    
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     misc.method = {'IE'};
%     formulation = {'BGU'};
%     misc.coreMethod = {'linear_FEM'};
%     degree = 1;
% 
%     collectIntoTasks
%     
%     misc.method = {'BA'};
%     formulation = {'VL2E'};
%     collectIntoTasks    
% end