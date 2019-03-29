

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';

method = {'IE'};
formulation = 'BGU';
coreMethod = 'C0IGA';

f = 1e3;             % Frequency
% f = 3.98e4;             % Frequency

M = 4; %3:6

alpha_s = 240*pi/180;
beta_s = 30*pi/180;
% alpha_s = pi;
% alpha_s = 0;
% beta_s = 0;
alpha = (0:0.5:360)*pi/180;
% alpha = 0;
% alpha = [0,90,180,270,360]*pi/180;
% beta = 30*pi/180;
beta = beta_s;
calculateVolumeError = 0;
calculateSurfaceError = 0;

calculateFarFieldPattern = 1;
plot3Dgeometry = 0;
degreeElev = 0;
parm = 1:2;

loopParameters = {'M','parm','method'};
collectIntoTasks
