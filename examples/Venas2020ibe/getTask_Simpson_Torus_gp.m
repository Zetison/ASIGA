scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
BC = 'NBC';
model = 'Torus';
method = {'BEM'};
formulation = {'CCBIE','CHBIE','CBM','CRCBIE1','CRCBIE2','CRCBIE3','GCBIE','GHBIE','GBM','GRCBIE1','GRCBIE2','GRCBIE3'};
% formulation = {'GBM'};
% formulation = {'GCBIE'};
% formulation = {'CCBIE','CBM'};
% formulation = {'CBM'};
% formulation = {'CCBIE'};
k = 0.2;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
M = [2,5]; 
% M = 1:2;
degree = [2,5];
% degree = [2,3,5];
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
% extraGPBEM = [0,1,2,3,4,5];
% extraGPBEM = 50;
extraGPBEM = [0,1,2,4,8,16,32];
extraGP = 0;
% extraGP = 4;
plot3Dgeometry = 0;  % Plot visualization of mesh and geometry in 3D
applyLoad = 'SimpsonTorus';
% beta = 0;
exteriorProblem = false;
agpBEM = 1:10;
% agpBEM = 1:3; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
% useNeumanProj = [1,0];
useNeumanProj = 0;
computeCondNumber = true;
plot3Dgeometry = 0;
quadMethodBEM = {'Simpson','New'};

solveForPtot = true;
loopParameters = {'extraGP','extraGPBEM','agpBEM','M','degree','method','formulation','quadMethodBEM'};
runTasksInParallel = 1;

collectIntoTasks


method = {'BA'};
quadMethodBEM = NaN;
formulation = {'SL2E'};
extraGPBEM = NaN;
extraGP = 0;
agpBEM = [agpBEM(1),agpBEM(end)]; % parameter for adaptiv Gauss point integration around singularities for BEM formulations

collectIntoTasks