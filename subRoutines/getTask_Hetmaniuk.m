scatteringCase = 'Sweep';
% BC = 'NNBC';
model = 'S1';  % Spherical shell

% coreMethod = {'IGA','hp_FEM'};
coreMethod = {'IGA'};
method = {'IE'};
formulation = {'BGU'};


% k = [9, 22.5, 36];
% k = 1;
% k = [9, 22.5, 36]/10;
k = [9, 22.5, 36]/5;
k = k(1):0.01:k(end);
c_f = 1500;
f = k*c_f/(2*pi);

alpha = 0;
beta = -pi/2;   
r_a = 1.2;

M = 5; % 5
N = 7; % 9

degree = 4;
% degree = 2:5;
parm = 1;
calculateSurfaceError = 1;
calculateVolumeError  = 0;
calculateFarFieldPattern = 0;
runTasksInParallel = true;
plot3Dgeometry = 0;
plot2Dgeometry = 0;
% initMeshFactXi = 3;
% initMeshFactZeta = 4;
useROM = false;
useDGP = false;
computeCondNumber = 0;

loopParameters = {'f','method','coreMethod','formulation','M','degree','parm'};
collectIntoTasks

method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_ROM = k(1):0.01:k(end);
% k_ROM = k(1):0.1:k(end);
% k_ROM = k(1):0.05:k(end);
% k_ROM = k(1):0.005:k(end);
% k_ROM = k(1):0.5:k(end);
% k_temp = k;
% k = k_ROM;
% k = (9:1:36)/10;
% k = 36;
runTasksInParallel = false;

method = {'IE'};
formulation = {'BGC'};
noVecsArr = [2,4,8,16,24,32,64];        % do not put noVecsArr in loopParameters (this is done automatically)
basisROMcell = {'DGP','Hermite','Bernstein','Pade','Taylor'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'Pade'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'Taylor'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'Bernstein'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'Hermite'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'DGP','Hermite'};  % do not put basisROMcell in loopParameters (this is done automatically)
useDGP = true;
% noVecsArr = [8,16];        % do not put noVecsArr in loopParameters (this is done automatically)
% noVecsArr = 1;
% noVecsArr = [8,32,64,128];
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
M = 5; % 7
N = 7; % 9
% M = 1; 
% N = 2;
useROM = true;
storeFullVarCol = true;
loopParameters = {'method','coreMethod','formulation','M','degree','parm'};
% collectIntoTasks

M = M-1; 
coreMethod = {'hp_FEM'};
% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
useROM = false;
coreMethod = {'hp_FEM'};
method = {'BA'};
formulation = {'VL2E'};
% collectIntoTasks

formulation = {'SL2E'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coreMethod = {'IGA'};
M = M+1; 
formulation = {'VL2E'};
% collectIntoTasks

formulation = {'SL2E'};
% collectIntoTasks