scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
model = 'M31';

parm = 1;
plot2Dgeometry = 0;
plot3Dgeometry = 0;

method = {'BEM'};
formulation = 'GCBIE';
solveForPtot = true;

BC = 'NNBC';
% BC = 'SSBC';
f = 1e2; 
M = 1:3; %1:5
alpha_s = 240*pi/180;
beta_s = 0;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
plotResultsInParaview = 0;

extraGP = 0; % extra quadrature points
degree = 2;
loopParameters = {'M','method'};
collectIntoTasks

% method = 'BA';
% method = {'IENSG'};
% method = {'BEM'};
% BC = {'SHBC', 'SSBC','NNBC'};
% for BC = {'SHBC', 'SSBC','NNBC'}
%     method = {'IE'};
% 
% %     coreMethod = {'C0_IGA','hp_FEM'};
% %     coreMethod = {'IGA'};
%     coreMethod = {'hp_FEM','h_FEM','C0_IGA','IGA'};
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
%     calculateSurfaceError = 1;
%     plot2Dgeometry = 0;
%     plot3Dgeometry = 0;
%     loopParameters = {'M','degree','method','coreMethod','formulation','BC'};
% 
%     N = 6;
%     
%     collectIntoTasks
%     
%     method = {'BA'};
%     formulation = {'VL2E'};
%     collectIntoTasks
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     method = {'IE'};
%     coreMethod = {'IGA'};
%     formulation = {'BGU'};
%     degree = 3:4;
% %     degree = 4;
%     collectIntoTasks
%     
%     method = {'BA'};
%     formulation = {'VL2E'};
%     collectIntoTasks    
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     method = {'IE'};
%     formulation = {'BGU'};
%     coreMethod = {'linear_FEM'};
%     degree = 1;
% 
%     collectIntoTasks
%     
%     method = {'BA'};
%     formulation = {'VL2E'};
%     collectIntoTasks    
% end