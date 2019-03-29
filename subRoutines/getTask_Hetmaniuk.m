scatteringCase = 'Sweep';
% BC = 'NNBC';
model = 'S1';  % Spherical shell

coreMethod = 'IGA';
% coreMethod = {'IGA'};
method = {'IE'};
formulation = {'BGU'};


k = 9:0.05:36;
% k = 36;
c_f = 1500;
f = k*c_f/(2*pi);

alpha_s = 0;
beta_s = -pi/2;  

alpha = 0;
beta = -pi/2;   
r_a = 1.2;

degree = 4;
M = 6; 
parm = 2;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
plot3Dgeometry = 0;
plot2Dgeometry = 0;
loopParameters = {'degree','M','parm','method','formulation'};

% collectIntoTasks


method = {'BA'};
formulation = {'SL2E'};
% collectIntoTasks

method = {'IE'};
formulation = {'BGC'};
noVecsArr = [8,16,24,32,64,128];        % do not put noVecsArr in loopParameters (this is done automatically)
basisROMcell = {'Bernstein','Taylor'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'Taylor'};  % do not put basisROMcell in loopParameters (this is done automatically)
% noVecsArr = [8,16,32];
% k_start = 9/10;
% k_end = 36/10;
% P = 3;
% n = P*noVecs;
% temp = cos(pi/(2*n));
% a = ((k_start+k_end)*temp-(k_end-k_start))/(2*temp);
% b = ((k_start+k_end)*temp+(k_end-k_start))/(2*temp);
% j = 1:n;
% k_arr = 1/2*(a+b)+1/2*(b-a)*cos((2*n-2*j+1)*pi/(2*n));
% k = k_arr(round(linspace(1,n,P)));
k = [9, 22.5, 36];
% k = [9, 22.5, 36]/10;
k_ROM = k(1):0.05:k(end);
% k_ROM = k(1):0.5:k(end);
c_f = 1500;
f = k*c_f/(2*pi);
M = 7; 
N = 9;
useROM = true;
storeFullVarCol = true;
loopParameters = {'degree','M','parm','method','formulation'};
collectIntoTasks

