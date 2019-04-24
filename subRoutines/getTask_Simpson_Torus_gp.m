scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
BC = 'NBC';
model = 'Torus';
method = {'BEM'};
formulation = {'CCBIE','CHBIE','CBM','GCBIE','GHBIE','GBM'};
% formulation = {'GBM'};
% formulation = {'GCBIE'};
formulation = {'CCBIE','CBM'};
% formulation = {'CBM'};
% formulation = {'CCBIE'};
k = 2;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
M = [2,5]; 
% M = 1:4;
degree = [2,5];
% degree = [2,3,5];
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
extraGPBEM = [0,5,10,15,20];
% extraGPBEM = 20;
extraGP = [0,1,2,3];
% extraGP = 4;
plot3Dgeometry = 0;  % Plot visualization of mesh and geometry in 3D
applyLoad = 'SimpsonTorus';
% beta = 0;
exteriorProblem = false;
agpBEM = 1:0.5:3; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
% agpBEM = 4; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
% useNeumanProj = [1,0];
useNeumanProj = 1;
computeCondNumber = true;
plot3Dgeometry        = 0;
 
loopParameters = {'M','degree','method','formulation','extraGP','extraGPBEM','agpBEM','colBEM_C0'};

% applyLoad = 'radialPulsation'; % with analytic solution for arbitrary geometries
collectIntoTasks
