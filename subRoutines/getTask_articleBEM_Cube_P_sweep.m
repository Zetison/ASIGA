scatteringCase = 'MS';

model = 'Cube_P';  % Spherical shell
BC = 'NBC';
coreMethod = 'IGA';
% coreMethod = {'IGA'};
method = {'BEM'};
formulation = {'CCBIE', 'CBM', 'CHBIE'};
% formulation = {'CCBIE'};


eigenValues = [1,2,3,4,5,6,8,9,10];
eigenValues = unique(sort([eigenValues, 3,6,9]));
eigenValues = pi*sqrt(eigenValues);%  Analytical eigenvalues of the interior Dirichlet/Neumann cube problem.

k = sort([eigenValues linspace(0.01,10,3000)]);
% k = sort([eigenValues linspace(0.01,10,100)]);
% k = 1;
f = k*1500/(2*pi);

alpha_s = 240*pi/180;
beta_s = 30*pi/180;  

alpha = 240*pi/180;
beta = 30*pi/180;   

applyLoad = 'radialPulsation'; % with analytic solution for arbitrary geometries
degree = 4;
M = 4;
plot3Dgeometry = 0;
calculateSurfaceError = 1;

runTasksInParallel = 1;
loopParameters = {'f','method','formulation'};

collectIntoTasks

method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks
