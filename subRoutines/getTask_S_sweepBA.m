scatteringCase = 'Sweep';

model = 'S5';  % Spherical shell

applyLoad = 'radialPulsation';
coreMethod = 'IGA';
method = 'BA';
formulation = 'SL2E';
calculateFarFieldPattern = false;
storeFullVarCol = true;

Nhat = 2^10;
f_c = 1500;
c_f = 1500;
T = 120/f_c;

% omega = linspace(1,Nhat/2-1,2^10+1)*2*pi/T; 
omega_min = 2*pi/T;
omega_max = (Nhat/2-1)*omega_min;
k_min = omega_min/c_f;
k_max = omega_max/c_f;
n = 10;
k = (k_min+k_max)/2 + (k_max-k_min)/2*cos((2*(n:-1:1).'-1)/(2*n)*pi);
omega = k*c_f;
% omega = 1*2*pi/T; 
f = omega/(2*pi);

alpha_s = pi;
beta_s = 0;  

alpha = alpha_s;
beta = beta_s;   

M = 6;
degreeElev = 4;
calculateSurfaceError = 1;
loopParameters = {'M'};
storeSolution = true;
parm = 1;

collectIntoTasks
