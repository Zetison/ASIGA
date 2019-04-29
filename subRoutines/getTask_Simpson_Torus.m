scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
BC = 'NBC';
model = 'Torus';
method = {'BEM'};
formulation = {'CCBIE','CHBIE','CBM','GCBIE','GHBIE','GBM'};
formulation = {'GCBIE','GBM'};
% formulation = {'GCBIE'};
% formulation = {'CCBIE'};
k = 2;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
M = 1:5;
% M = 3; 
degree = 2:5;
% degree = 2;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
% extraGPBEM = [0,3];
% extraGP = [0,3];
extraGPBEM = [4,8,32];
% extraGPBEM = 32;
% extraGP = 2;
extraGP = 1;
plot3Dgeometry = 0;  % Plot visualization of mesh and geometry in 3D
applyLoad = 'SimpsonTorus';
% beta = 0;
exteriorProblem = false;
% agpBEM = [1,2,4,8]; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
% agpBEM = 2; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
agpBEM = [2,4,8]; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
% useNeumanProj = [1,0];
useNeumanProj = 0;
computeCondNumber = 0;
runTasksInParallel = false;

loopParameters = {'M','degree','method','formulation','extraGP','extraGPBEM','agpBEM','useNeumanProj'};

% applyLoad = 'radialPulsation'; % with analytic solution for arbitrary geometries
collectIntoTasks


useNeumanProj = [1,0];
formulation = {'CCBIE','CBM'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
useNeumanProj = 0;
% degree = 2;
extraGPBEM = 4;
extraGP = 4;
agpBEM = 2; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
M = 1:7;
method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks