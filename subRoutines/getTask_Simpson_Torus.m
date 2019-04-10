scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
BC = 'NBC';
model = 'Torus';
method = {'BEM'};
formulation = {'GCBIE','CBM','CCBIE','GBM'};
formulation = {'GBM'};
k = 2;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
M = 1:5;
% M = 1:2; 
degree = 2:4;
degree = 2;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
% extraGPBEM = [0,4];
% extraGP = [0,4];
extraGPBEM = 5;
extraGP = 5;
plot3Dgeometry = 0;  % Plot visualization of mesh and geometry in 3D
applyLoad = 'SimpsonTorus';
% beta = 0;
exteriorProblem = false;
agpBEM = [1,2,4,8]; % parameter for adaptiv Gauss point integration around singularities for BEM formulations

loopParameters = {'method','formulation','M','degree','extraGP','extraGPBEM','agpBEM'};

% applyLoad = 'radialPulsation'; % with analytic solution for arbitrary geometries
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
% degree = 2;
agpBEM = 2; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks