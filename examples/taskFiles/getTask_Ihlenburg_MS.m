function studies = getTask_Ihlenburg_MS()
% This study is based on Ihlenburg1998fea (Figure 7)
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271
% Experimental data at scale 1:50
counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
scatteringCase = 'Sweep';
BC = 'SSBC';
model = 'IMS';  % Spherical shell

coreMethod = {'IGA','hp_FEM'};
coreMethod = {'IGA'};
method = {'IE'};
applyLoad = 'pointSource';
formulation = {'PGC'};

varCol{1} = struct('media', 'fluid', ...
                   't', [0.1143,0.04064], ...
                   'R1', 4.6863, ...
                   'R2', 4.6863, ...
                   'L', 67.9958, ...
                   'c_f', 1482, ...
                   'rho', 1000);
varCol{2} = struct('media', 'solid', ...
                   'E', 2.0e11, ...
                   'nu', 0.29, ...
                   'rho', 7908.5);

varCol{1}.meshFile = 'createNURBSmesh_M3';
warning('off','NURBS:weights')

k = linspace(2.5, 20, 5)/varCol{1}.R1;
% k = [9, 22.5, 36]/5;

alpha = 0;
beta = 0;   
sdfmsc = 3*(varCol{1}.R1+varCol{1}.L+varCol{1}.R2)/1.36;
r_s = sdfmsc - varCol{1}.L/2;

M = 1:5; % 5
M = 3; % 5
N = 1; % 9

degree = 4;
% degree = 2:5;
parm = 1;
calculateSurfaceError = 0;
calculateVolumeError  = 0;
calculateFarFieldPattern = 0;
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 1;
computeCondNumber = 0;

postPlot(1).xname        	= 'k_ROM';
postPlot(1).yname        	= 'TS';
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

k_ROM = k(1):2.5:k(end);
c_f = 1500;
f = k*c_f/(2*pi);

noVecsArr = [4,8,16];        % do not put noVecsArr in loopParameters (this is done automatically)
% basisROMcell = {'Bernstein','Hermite','Pade','Taylor','DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'Pade'};  % do not put basisROMcell in loopParameters (this is done automatically)
% basisROMcell = {'Taylor'};  % do not put basisROMcell in loopParameters (this is done automatically)
basisROMcell = {'DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
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
N = 13; % 9
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
p_ie = 5;
s_ie = 2;
IElocSup = 0;        % Toggle usage of radial shape functions in IE with local support

storeFullVarCol = false;
loopParameters = {'method','IElocSup'};
% collectIntoTasks

N = 40; % 9
IElocSup = 1; 
% collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
useROM = false;
postPlot(1).xname = 'k';
postPlot(1).xLoopName = 'k';
loopParameters = {'method'};
k = k_ROM;
c_f = 1500;
f = k*c_f/(2*pi);
coreMethod = {'IGA'};
method = {'BA'};
formulation = {'SL2E'};
% collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = {'IE'};
formulation = {'PGC'};
N = 13; % 9
IElocSup = 0;
loopParameters = {'method','IElocSup'};
collectIntoTasks
p_ie = 5;
N = p_ie*5; % 9
IElocSup = 1;
% collectIntoTasks
