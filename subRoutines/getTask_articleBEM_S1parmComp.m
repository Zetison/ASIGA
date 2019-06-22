scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
model = 'S1';
method = {'IE'};
formulation = {'BGU'};
k = 1;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
N = 6;
M = 1:6;
% M = 1;
parm = [1,2];
% parm = 2;
degree = 4;

alpha_s = 240*pi/180;
beta_s = 30*pi/180;

alpha = (0:0.5:360)*pi/180;
beta = 30*pi/180;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;

loopParameters = {'M','parm','method','formulation'};

% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
method = {'BA'};
formulation = {'SL2Etot','SL2E'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = {'BEM'};
formulation = {'CCBIE', 'CBM', 'CHBIE', 'GCBIE', 'GBM', 'GHBIE'};
useNeumannProj = false;

% collectIntoTasks
