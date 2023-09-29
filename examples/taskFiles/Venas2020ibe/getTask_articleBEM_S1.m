misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
misc.model = 'S1';
misc.method = {'IE'};
formulation = {'BGU'};
k = 1;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
N = 6;
M = 1:6;
parm = [1,2];
degree = 4;
alpha_s = 240*pi/180;
beta_s = 30*pi/180;
% alpha_s = 0;
% beta_s = 0;
alpha = alpha_s;
beta = beta_s;
err.calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
prePlot.plot3Dgeometry        = 0; 
extraGPBEM = 0;
extraGP = 0;

% beta = 0;

solveForPtot = true;
loopParameters = {'M','degree','misc.method','formulation','parm'};

% misc.applyLoad = 'pointPulsation'; % with analytic solution for arbitrary geometries
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
misc.method = {'BA'};
formulation = {'SL2E'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

misc.method = {'BEM'};
M = 1;
degree = 2;
parm = 1;
extraGPBEM = 0;
extraGP = 0;
% formulation = {'CRCBIE1', 'CRCBIE2', 'CRCBIE3', 'GRCBIE1', 'GRCBIE2', 'GRCBIE3', 'CCBIE', 'GCBIE', 'CHBIE', 'GHBIE', 'CBM', 'GBM'};
% formulation = {'CBM','GBM'};
formulation = {'CCBIE'};
formulation = {'CHBIE'};
collectIntoTasks
