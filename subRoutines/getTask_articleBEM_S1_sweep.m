scatteringCase = 'MS';

model = 'S1';  % Spherical shell

coreMethod = 'IGA';
% coreMethod = {'IGA'};
method = {'BEM'};
formulation = {'CCBIE', 'CBM', 'CHBIE', 'CCBIEC'};
% formulation = {'CCBIE'};
solveForPtot = true;


% Eigenfrequencies from Zheng2015itb
eigenValues = [  pi  % Analytical eigenvalues of the interior Dirichlet sphere example.
               2*pi
               3*pi
               4.493409457909065
               7.725251836937707
               5.763459196894550
               9.095011330476353
               6.987932000500519
               8.182561452571243
               9.355812111042747
               4.493409457909064 % Analytical eigenvalues of the interior Neumann sphere example.
               7.725251836937708
               2.081575977818101
               5.940369990572713
               9.205840142936665
               3.342093657365695
               7.289932304093351
               4.514099647032276
               8.583754956365768
               5.646703620436797
               9.840446043040137               
               6.756456330204129
               7.851077679474405
               8.934838878352839
               ]';
if false 
    k = eigenValues;
else
    noPts = 1000;
%     noPts = 100;
    delta = 10/noPts*3;
    k = linspace(0.01,10,noPts);
    k = [k, linspace(0.01,0.05,round(noPts/10))];
    for i = 1:numel(eigenValues)
        k = [k, eigenValues(i)+linspace(-delta/2,delta/2,round(noPts/10))];
    end
end
k = sort(unique(k));
% k = sort([eigenValues linspace(0.005,10,3000)]);
% k = sort([eigenValues linspace(0.005,10,100)]);
% k = 1;
f = k*1500/(2*pi);

alpha_s = 240*pi/180;
beta_s = 30*pi/180;  

alpha = 240*pi/180;
beta = 30*pi/180;   

degree = 4;
M = 4;
parm = 2;
calculateSurfaceError = 1;

runTasksInParallel = 1;
loopParameters = {'f','method','formulation'};

collectIntoTasks

method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks

method = {'IE'};
formulation = {'BGU'};
N = 5;
solveForPtot = false;
% collectIntoTasks
