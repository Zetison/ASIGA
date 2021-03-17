misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
misc.model = 'S1';
k = 1;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
N = 6;
M = 1:6;
M = 1:2;
parm = [1,2];
parm = 2;
degree = 4;

alpha_s = 240*pi/180;
beta_s = 30*pi/180;

alpha = (0:0.5:360)*pi/180;
beta = 30*pi/180;
err.calculateSurfaceError = 1;
calculateSurfEnrgErr = true;
calculateFarFieldPattern = 1;
solveForPtot = 1;
% runTasksInParallel = true;

solveForPtot = true;
loopParameters = {'M','parm','misc.method','formulation'};

misc.method = {'BEM'};
formulation = {'CCBIE', 'CBM', 'CHBIE', 'GCBIE', 'GBM', 'GHBIE'};
formulation = {'CCBIE', 'GCBIE'};
% formulation = {'CCBIE'};
useNeumannProj = false;

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
misc.method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks
