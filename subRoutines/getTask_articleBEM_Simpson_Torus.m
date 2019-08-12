scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
BC = 'NBC';
model = 'Torus';
method = {'BEM'};
formulation = {'GCBIE','GBM','CCBIE','CBM'};
k = 2;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
M = 1:5; %1:5
degree = [2,4];
calculateSurfaceError = 1;
calculateSurfEnrgErr = true;
calculateFarFieldPattern = 0;

plot3Dgeometry = 0;  % Plot visualization of mesh and geometry in 3D
applyLoad = 'SimpsonTorus';
exteriorProblem = false;
useNeumanProj = 0;
computeCondNumber = 0;
runTasksInParallel = false;

solveForPtot = false;
loopParameters = {'M','degree','method','formulation','extraGP','extraGPBEM','agpBEM','useNeumanProj'};

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
useNeumanProj = 0;
extraGP = 0;
M = 1:7;
method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks