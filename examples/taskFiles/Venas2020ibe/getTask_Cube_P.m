misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
BC = 'NBC';
misc.model = 'Cube_P';
misc.method = {'BEM'};
formulation = {'CCBIE','CBM'};
k = 2;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
M = 1:5;
% M = 1; 
% degree = 2:5;
degree = [1,2,4];
% degree = 2;
err.calculateSurfaceError = 1;
calculateFarFieldPattern = 0;

prePlot.plot3Dgeometry = 0;  % Plot visualization of mesh and geometry in 3D
misc.applyLoad = 'pointPulsation';
exteriorProblem = true;
computeCondNumber = 0;
runTasksInParallel = 0;
useNeumanProj = [1,0];

solveForPtot = false;
loopParameters = {'M','degree','misc.method','formulation','extraGP','extraGPBEM','agpBEM','useNeumanProj','colBEM_C0'};
colBEM_C0 = [0,1/2];
colBEM_C0 = 0;
useNeumanProj = 0;
% degree = 2;
% M = 1:5;
calculateSurfEnrgErr = false;

collectIntoTasks

formulation = {'GCBIE','GBM'};
useNeumanProj = 0;
colBEM_C0 = NaN;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
% degree = 2;
extraGP = 1;
M = 1:7;
misc.method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks
