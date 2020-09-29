function studies = getTask_Hetmaniuk()
% This study is based on Hetmaniuk2012raa (Figure 7)
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271

counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
scatteringCase = 'Sweep';
% BC = 'NNBC';
model = 'S1';  % Spherical shell

coreMethod = {'IGA','hp_FEM'};
coreMethod = {'IGA'};
method = {'IE'};
formulation = {'BGC'};

varCol = setS1Parameters('double',1);
varCol{1}.meshFile = 'createNURBSmesh_EL';

% k = [9, 22.5, 36];
k = [9, 22.5, 36]/5;

alpha = 0;
beta = -pi/2;   
r_a = 1.2;

M = 1:5; % 5
N = 1; % 9

degree = 4;
% degree = 2:5;
parm = 2;
calculateSurfaceError = 1;
calculateVolumeError  = 0;
calculateFarFieldPattern = 0;
plot3Dgeometry = 0;
plot2Dgeometry = 0;
% initMeshFactXi = 3;
% initMeshFactZeta = 4;
computeCondNumber = 0;

postPlot(1).xname        	= 'k_ROM';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'semilogy';
postPlot(1).lineStyle    	= '-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).xLoopName     	= 'k_ROM';
postPlot(1).legendEntries 	= {};
postPlot(1).subFolderName 	= '';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;

k_ROM = k(1):0.05:k(end);
c_f = 1500;
f = k*c_f/(2*pi);

noVecsArr = [8,16,24,32,64];        % do not put noVecsArr in loopParameters (this is done automatically)
basisROMcell = {'Bernstein','Hermite','Pade','Taylor','DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'Pade'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'Taylor'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'Hermite'};  % do not put basisROMcell in loopParameters (this is done automatically)
% noVecsArr = 8;
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
M = 4; % 5
N = 7; % 9
% M = 1; 
% N = 2;
useROM = true;
% 
% useROM = false;
% postPlot(1).xname = 'k';
% k = k_ROM;
% c_f = 1500;
% f = k*c_f/(2*pi);

storeFullVarCol = false;
loopParameters = {'method','coreMethod','formulation','M','degree','parm'};
collectIntoTasks

M = M-1; 
coreMethod = {'hp_FEM'};
% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
useROM = false;
postPlot(1).xname = 'k';
k = k_ROM;
c_f = 1500;
f = k*c_f/(2*pi);
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