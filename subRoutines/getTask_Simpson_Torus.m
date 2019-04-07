scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
BC = 'NBC';
model = 'Torus';
method = {'BEM'};
formulation = {'CBM','CCBIE','GBM','GCBIE'};
k = 2;
% f = 1e2;             % Frequency
f = k*1500/(2*pi);
M = 1:5; % M = 1 -> ndof = 62?
degree = 2:4;
% degree = 2:4;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
extraGPBEM = 4;
extraGP = 4;
plot3Dgeometry        = 0;  % Plot visualization of mesh and geometry in 3D
applyLoad = 'SimpsonTorus';
% beta = 0;
exteriorProblem = false;

loopParameters = {'M','degree','method','formulation'};

% applyLoad = 'radialPulsation'; % with analytic solution for arbitrary geometries
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
% degree = 2;
method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks