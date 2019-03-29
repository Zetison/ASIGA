scatteringCase = 'Sweep';

model = 'S1';  % Spherical shell

coreMethod = 'IGA';
% coreMethod = {'IGA'};
method = {'BEM'};
if strcmp(method, 'BEM')
    formulation = {'CCBIE', 'CBM', 'CHBIE'};
end


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
k = sort([eigenValues linspace(0.005,10,3000)]);
f = k*1500/(2*pi);

alpha_s = 240*pi/180;
beta_s = 30*pi/180;  
alpha_s = 0;
beta_s = 0;  

alpha = alpha_s;
beta = beta_s;   

M = 1:3;
calculateSurfaceError = 1;
loopParameters = {'M','method','formulation'};

% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
f = linspace(1e3,3e5,3000); 
% f = [1e4 2e4]; 
method = {'MFS'};
coreMethod = 'IGA';
formulation = '';
M = 5;
degreeElev = 0;
calculateSurfaceError = 0;
computeCondNumber = false;
loopParameters = {'M','method'};
% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT simulation
method = {'RT'};
formulation = '';
M = 3;
plot3Dgeometry = 0;
degreeElev = 0;
calculateSurfaceError = 0;
computeCondNumber = false;
plotFarField = 1;
applyLoad = 'planeWave';
parm = 3:6;
r = 2;

loopParameters = {'parm','M','method'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KDT simulation
method = {'KDT'};
formulation = '';
M = 3;
plot3Dgeometry = 0;
degreeElev = 0;
calculateSurfaceError = 0;
computeCondNumber = false;
loopParameters = {'M','method'};
% collectIntoTasks

