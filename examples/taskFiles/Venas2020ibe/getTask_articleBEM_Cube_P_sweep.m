misc.scatteringCase = 'MS';

misc.model = 'Cube_P';  % Spherical shell
BC = 'NBC';
misc.coreMethod = 'IGA';
% misc.coreMethod = {'IGA'};
misc.method = {'BEM'};
formulation = {'CCBIE', 'CBM', 'CHBIE', 'CCBIEC'};
% formulation = {'CCBIEC'};


eigenValues = [];
eigenValues = [eigenValues, 1,2,3,4,5,6,8,9,10];
eigenValues = [eigenValues, 3,6,9];
% eigenValues = [eigenValues, 9];
eigenValues = unique(sort(eigenValues));
eigenValues = pi*sqrt(eigenValues);%  Analytical eigenvalues of the interior Dirichlet/Neumann cube problem.
noPts = 1000;
%     noPts = 100;
delta = 10/noPts*3;
k = linspace(0.01,10,noPts);
k = [k, linspace(0.01,delta/2,round(noPts/10))];
for i = 1:numel(eigenValues)
    k = [k, eigenValues(i)+linspace(-delta/2,delta/2,round(noPts/10))];
end
k = sort(unique(k));
% k = pi*sqrt(3);
% k = pi*sqrt(6);
misc.omega = k*1500;

ffp.alpha = 240*pi/180;
ffp.beta = 30*pi/180;   

misc.applyLoad = 'pointPulsation'; % with analytic solution for arbitrary geometries
msh.degree = 4;
msh.M = 4;
prePlot.plot3Dgeometry = 0;
err.calculateSurfaceError = 1;

runTasksInParallel = 1;
misc.solveForPtot = false;
loopParameters = {'f','misc.method','formulation'};

collectIntoTasks

misc.method = {'BA'};
misc.formulation = {'SL2E'};
collectIntoTasks
