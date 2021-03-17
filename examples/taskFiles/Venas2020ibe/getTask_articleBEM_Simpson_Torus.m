misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
BC = 'NBC';
misc.model = 'Torus';
misc.method = {'BEM'};
formulation = {'GCBIE','GBM','CCBIE','CBM'};
k = 2;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
M = 1:6; %1:5
degree = 2:4;
err.calculateSurfaceError = 1;
calculateSurfEnrgErr = true;
calculateFarFieldPattern = 0;

prePlot.plot3Dgeometry = 0;  % Plot visualization of mesh and geometry in 3D
misc.applyLoad = 'SimpsonTorus';
exteriorProblem = false;
useNeumanProj = 0;
computeCondNumber = 0;
runTasksInParallel = false;

solveForPtot = false;
loopParameters = {'M','degree','misc.method','formulation','extraGP','extraGPBEM','agpBEM','useNeumanProj'};

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
useNeumanProj = 0;
extraGP = 0;
M = 1:6;
misc.method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks