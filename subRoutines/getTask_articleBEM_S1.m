scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
model = 'S1';
method = {'IE'};
formulation = {'BGU'};
k = 1;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
N = 6;
M = 1:7;
parm = 1;
degree = 2;
alpha_s = 240*pi/180;
beta_s = 30*pi/180;
% alpha_s = 0;
% beta_s = 0;
alpha = (0:0.5:360)*pi/180;
beta = 30*pi/180;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
extraGPBEM = 0;
extraGP = 0;

% beta = 0;

loopParameters = {'M','degree','method','formulation','parm'};

% applyLoad = 'radialPulsation'; % with analytic solution for arbitrary geometries
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method = {'BEM'};
M = 1:6;
% coreMethod = {'IGA', 'C0_IGA', 'hp_FEM','h_FEM'};
% coreMethod = {'linear_FEM'};
% formulation = {'CCBIE', 'CBM', 'CHBIE', 'GCBIE', 'GBM', 'GHBIE'};
formulation = {'CCBIE', 'CBM', 'CHBIE', 'GCBIE', 'GBM', 'GHBIE','GRCBIE','CRCBIE'};
% formulation = {'CCBIE', 'CBM', 'CHBIE'};
% formulation = {'CRCBIE'};
% formulation = {'GRCBIE','CRCBIE','GCBIE','CCBIE'};
% formulation = {'CCBIE'};
% formulation = {'GCBIE'};
% formulation = {'CCBIE', 'GCBIE', 'GBM'};

collectIntoTasks
