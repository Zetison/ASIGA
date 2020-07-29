scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
BC = 'NBC';
model = 'Cube';
method = {'BEM'};
formulation = {'CCBIE','CHBIE','CBM','GCBIE','GHBIE','GBM'};
% formulation = {'GCBIE','GBM'};
% formulation = {'GCBIE'};
% formulation = {'CCBIE','CHBIE','CBM','CRCBIE1','CRCBIE2','CRCBIE3'};
% formulation = {'CRCBIE1','CRCBIE2','CRCBIE3'};
% formulation = {'CCBIE','CRCBIE1','CRCBIE2','CRCBIE3'};
k = 2;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
M = 1:3;
% M = 3; 
% degree = 2:5;
degree = [1,5];
degree = 2;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
% extraGPBEM = [0,3];
% extraGP = [0,3];
% extraGPBEM = [4,8,32];
extraGPBEM = 50;
% extraGP = 2;
extraGP = 0;
prePlot.plot3Dgeometry = 0;  % Plot visualization of mesh and geometry in 3D
applyLoad = 'SimpsonTorus';
% beta = 0;
exteriorProblem = false;
% agpBEM = [1,2,4,8]; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
agpBEM = 6; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
% agpBEM = [2,4,8]; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
% useNeumanProj = [1,0];
computeCondNumber = 0;
runTasksInParallel = false;
useNeumanProj = [1,0];

solveForPtot = false;
loopParameters = {'M','degree','method','formulation','extraGP','extraGPBEM','agpBEM','useNeumanProj'};
colBEM_C0 = 0;
quadMethodBEM = 'New';

% applyLoad = 'pointPulsation'; % with analytic solution for arbitrary geometries
% collectIntoTasks

formulation = {'GCBIE','GHBIE','GBM','GRCBIE1','GRCBIE2','GRCBIE3'};
useNeumanProj = 0;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
useNeumanProj = 0;
% degree = 2;
extraGPBEM = NaN;
extraGP = 1;
agpBEM = NaN; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
% M = 1:7;
method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks
