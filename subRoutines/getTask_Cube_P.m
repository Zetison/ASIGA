scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
BC = 'NBC';
model = 'Cube_P';
method = {'BEM'};
% formulation = {'CCBIE','CHBIE','CBM','GCBIE','GHBIE','GBM'};
% formulation = {'GCBIE','GBM'};
% formulation = {'GCBIE'};
% formulation = {'CCBIE','CHBIE','CBM','CRCBIE1','CRCBIE2','CRCBIE3'};
% formulation = {'CRCBIE1','CRCBIE2','CRCBIE3'};
formulation = {'CCBIE','CBM'};
k = 2;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
M = 1:5;
% M = 1; 
% degree = 2:5;
degree = 1:4;
% degree = 2;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
% extraGPBEM = [0,3];
% extraGP = [0,3];
% extraGPBEM = [4,8,32];
extraGPBEM = 16;
% extraGP = 2;
extraGP = 0;
plot3Dgeometry = 0;  % Plot visualization of mesh and geometry in 3D
applyLoad = 'radialPulsation';
% beta = 0;
exteriorProblem = true;
% agpBEM = [1,2,4,8]; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
agpBEM = 1.4; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
% agpBEM = [2,4,8]; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
% useNeumanProj = [1,0];
computeCondNumber = 0;
runTasksInParallel = 0;
useNeumanProj = [1,0];

loopParameters = {'M','degree','method','formulation','extraGP','extraGPBEM','agpBEM','useNeumanProj','colBEM_C0'};
colBEM_C0 = [0,1/2];
quadMethodBEM = 'Adaptive';

collectIntoTasks

formulation = {'GCBIE','GBM'};
useNeumanProj = 0;
colBEM_C0 = NaN;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
% degree = 2;
extraGPBEM = NaN;
extraGP = 1;
agpBEM = NaN; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
M = 1:7;
method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks