

scatteringCase = {'BI'};

model = {'S3'};
coreMethod = {'XI'};
% coreMethod = {'IGA'};

method = {'BEM'};
if strcmp(method, 'BEM')
    formulation = {'CBM'};
%     formulation = {'CCBIE', 'CBM'};
end

f_arr = 1e3;             % Frequency

M = 2;
degreeElevArr = 0;
alpha_s = 180*pi/180;
beta_s = 0*pi/180;   
% alpha_s = 240*pi/180;
% beta_s = 30*pi/180;          
% plot3Dgeometry = true;
calculateSurfaceError = 1;	% Only for spherical shell and if scatteringCase == 'Bi'
solveForPtot = true;

collectIntoTasks