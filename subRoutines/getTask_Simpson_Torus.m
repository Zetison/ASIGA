scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
BC = 'NBC';
model = 'Torus';
method = {'BEM'};
formulation = {'GCBIE','CBM','CCBIE','GBM'};
formulation = {'GBM'};
formulation = {'GCBIE'};
k = 2;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
M = 1:5;
% M = 5; 
degree = 2:4;
degree = 4;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
% extraGPBEM = [0,4];
% extraGP = [0,4];
extraGPBEM = 4;
extraGP = 4;
plot3Dgeometry = 0;  % Plot visualization of mesh and geometry in 3D
applyLoad = 'SimpsonTorus';
% beta = 0;
exteriorProblem = false;
% agpBEM = [1,2,4,8]; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
agpBEM = 2; % parameter for adaptiv Gauss point integration around singularities for BEM formulations

loopParameters = {'method','formulation','M','degree','extraGP','extraGPBEM','agpBEM'};

% applyLoad = 'radialPulsation'; % with analytic solution for arbitrary geometries
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
% degree = 2;
agpBEM = 2; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
M = 1:7;
degree = [2,3,4];
method = {'BA'};
formulation = {'SL2E'};
% collectIntoTasks