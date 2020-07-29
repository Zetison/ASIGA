scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
BC = 'NBC';
model = 'Torus';
method = {'BEM'};
formulation = {'GCBIE','GBM','CCBIE','CBM'};
k = 2;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
M = 1:5;
degree = 2:5;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;

prePlot.plot3Dgeometry = 0;  % Plot visualization of mesh and geometry in 3D
applyLoad = 'SimpsonTorus';
exteriorProblem = false;
useNeumanProj = 0;
computeCondNumber = 0;
runTasksInParallel = false;

loopParameters = {'M','degree','method','formulation','extraGP','extraGPBEM','agpBEM','useNeumanProj'};

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
useNeumanProj = 0;
extraGPBEM = NaN;
extraGP = 0;
agpBEM = NaN; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
M = 1:7;
method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks