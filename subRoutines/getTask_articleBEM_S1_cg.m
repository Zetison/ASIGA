scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';


% method = {'BA'};
% formulation = {'SL2E'};
method = {'BEM'};
formulation = {'GCBIE'};
coreMethod = {'IGA'};

c_f = 1500;
k = 1;
omega = k*c_f;
f = omega/(2*pi); 
% f = 1000;
% omega = 2*pi*f;
% k = omega/c_f;

M = 2:3;

alpha_s = 240*pi/180;
beta_s = 30*pi/180;

parm = 2;
degree = 4:5;
plotResultsInParaview = 0;
calculateFarFieldPattern = 0;
calculateSurfaceError = 1;
plot2Dgeometry = 0;
plot3Dgeometry = 0;

% extraGP = 1;
% extraGPBEM = 1;
plotResidualError = true;
loopParameters = {'M','formulation','method','degree','coreMethod'};
collectIntoTasks

Mend = 6;
degree = 4:5;
M = 1:Mend;
plotResidualError = false;
formulation = {'CCBIE','CHBIE','CBM'};
loopParameters = {'M','formulation','method','degree','coreMethod','colMethod'};
colMethod = {'CG','Grev','GL'};
collectIntoTasks


method = {'BA'};
formulation = {'SL2E'};
colMethod = NaN;
collectIntoTasks

