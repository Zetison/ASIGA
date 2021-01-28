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
% formulation = {'PGC'};
formulation = {'BGC'};
noDomains = 1;
varCol = setHetmaniukParameters(1);
varCol{1}.meshFile = 'createNURBSmesh_EL';
varCol{1}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, max(2^(M-1)/8-1,0)];
if noDomains > 1
    varCol{2}.refinement = @(M,t,t_fluid) [2^(M-1)-1, 2^(M-1)-1, max(round(t/t_fluid)*2^(M-1),0)];
end

k = linspace(9, 36, 3);
k = linspace(9, 36, 3)/5;

alpha = 0;
beta = -pi/2;   
r_a = 1.2;

parm = 1;
calculateSurfaceError = 1;
calculateVolumeError  = 0;
calculateFarFieldPattern = 1;
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;
prePlot.colorFun = @(v) abs(norm2(v)-1);
prePlot.resolution = [200,200,0];
computeCondNumber = 0;

postPlot(1).xname        	= 'f_ROM';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'semilogy';
postPlot(1).lineStyle    	= '-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).legendEntries 	= {};
postPlot(1).subFolderName 	= '';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;

postPlot(2) = postPlot(1);
postPlot(2).yname        	= 'TS';
postPlot(2).axisType      	= 'plot';

k_ROM = k(1):0.05:k(end);
k_ROM = k(1):0.1:k(end);
c_f = varCol{1}.c_f;
omega_ROM = k_ROM*c_f;
omega = k*c_f;
f = k*c_f/(2*pi);

% noVecsArr = [4,8,16];        % do not put noVecsArr in loopParameters (this is done automatically)
% basisROMcell = {'Bernstein','Hermite','Pade','Taylor','DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
basisROMcell = {'Pade','Taylor','DGP','Hermite','Bernstein'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'Taylor'};  % do not put basisROMcell in loopParameters (this is done automatically)
basisROMcell = {'DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'Hermite'};  % do not put basisROMcell in loopParameters (this is done automatically)
noVecsArr = 64;
% noVecsArr = [2,4,8,16,32,64];
% noVecsArr = 8;
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
degree = 2;
% degree = 2:5;
M = 3; % 5
N = 2; % 9
% M = 1; 
% N = 2;
useROM = true;
% 
% useROM = false;
% postPlot(1).xname = 'k';
% k = k_ROM;
% c_f = 1500;
% f = k*c_f/(2*pi);
Xi = [0,0,0,1,1,2,2,3,3,3]/3;
refineThetaOnly = true;
% refineThetaOnly = 0;
p_ie = 5;
s_ie = 2;
IElocSup = 0;        % Toggle usage of radial shape functions in IE with local support

storeFullVarCol = false;
loopParameters = {'M','method'};
% collectIntoTasks

N = 40; % 9
IElocSup = 0; 
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
useROM = false;
postPlot(1).xname = 'f';
postPlot(2).xname = 'f';
omega = omega_ROM;
f = omega/(2*pi);
coreMethod = {'IGA'};
method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = {'IE'};
formulation = {'PGC'};
N = 13; % 9
IElocSup = 0;
loopParameters = {'method','IElocSup'};
% collectIntoTasks
p_ie = 5;
N = p_ie*5; % 9
IElocSup = 1;
% collectIntoTasks
