

scatteringCase = 'Sweep'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1'; % BeTSSi model 5A

method = {'IE'};
% method = {'IENSG'};
formulation = {'BGU'};

f = linspace(1e2,2e2,2); %[1e2 5e2 1e3];             % Frequency
f = 1e3; %[1e2 5e2 1e3];             % Frequency
omega = 2*pi*f;
k = omega/1500;
N = 3;
parm = 1;

M = 1:4;
degree = 2;
loopParameters = {'M','degree','N','formulation','method'};
alpha = 240*pi/180;

alpha_s = 240*pi/180;
beta_s = 0*pi/180;        
% beta_s = 30*pi/180;
plot3Dgeometry = 0;
calculateFarFieldPattern = 0;
calculateSurfaceError = 1;

plot2Dgeometry = 0;  % Plot cross section of mesh and geometry
collectIntoTasks


formulation = {'BGC'};
N = 9;
collectIntoTasks
