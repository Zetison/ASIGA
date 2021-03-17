misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
BC = 'SHBC';
misc.model = 'Cube';
misc.method = {'BEM'};
formulation = {'CCBIE','CHBIE','CBM','GCBIE','GHBIE','GBM'};
% formulation = {'GCBIE','GBM'};
% formulation = {'GCBIE'};
% formulation = {'CCBIE','CHBIE','CBM','CRCBIE1','CRCBIE2','CRCBIE3'};
formulation = {'CCBIE','CHBIE','CBM','CRCBIE1','CRCBIE2','CRCBIE3','GCBIE','GHBIE','GBM','GRCBIE1','GRCBIE2','GRCBIE3'};
% formulation = {'CRCBIE1'};
formulation = {'CCBIE'};
k = 2;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
M = 1:3;
M = 1; 
% degree = 2:5;
degree = 2;
err.calculateSurfaceError = strcmp(BC,'NBC');
calculateFarFieldPattern = strcmp(BC,'SHBC');
alpha = (0:0.5:360)*pi/180;
beta = 0;
alpha_s = 240*pi/180;
beta_s = 30*pi/180;
% extraGPBEM = [0,3];
% extraGP = [0,3];
% extraGPBEM = [4,8,32];
extraGPBEM = 32;
% extraGP = 2;
extraGP = 0;
prePlot.plot3Dgeometry = 1;  % Plot visualization of mesh and geometry in 3D
if strcmp(BC,'SHBC')
    misc.applyLoad = 'planeWave';
else
    misc.applyLoad = 'pointPulsation';
end
% beta = 0;
exteriorProblem = true;
% agpBEM = [1,2,4,8]; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
agpBEM = 5; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
% agpBEM = [2,4,8]; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
% useNeumanProj = [1,0];
computeCondNumber = 0;
runTasksInParallel = false;

solveForPtot = true;
loopParameters = {'M','degree','misc.method','formulation','extraGP','extraGPBEM','agpBEM'};
% colBEM_C0 = 0;
colBEM_C0 = 1/2;
quadMethodBEM = 'New';

% misc.applyLoad = 'pointPulsation'; % with analytic solution for arbitrary geometries
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
useNeumanProj = 0;
% degree = 2;
extraGPBEM = NaN;
extraGP = 1;
agpBEM = NaN; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
% M = 1:7;
misc.method = {'BA'};
formulation = {'SL2E'};
% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
misc.method = {'MFS'};
formulation = 'PS'; % PS = point solution, SS = spherical solution
M = 1;
% parm = linspace(0.1,0.4,8);
% parm = 0.5/sqrt(k);
computeCondNumber = false;
solveForPtot = false;
% delta = 0.1;
delta = [0.1,0.2,0.4,0.6,1]/10;
extraGP = [1,2,4,8,16,32];
loopParameters = {'delta','M','extraGP','misc.method'};

% collectIntoTasks
