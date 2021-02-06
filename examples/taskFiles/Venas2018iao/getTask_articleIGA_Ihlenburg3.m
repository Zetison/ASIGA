function studies = getTask_articleIGA_Ihlenburg3()
% This study correspond to Figure 18 in Venas2018iao
% Venas2018iao is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)
counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'IL';

method = {'IE'};
formulation = 'BGU';

BC = 'NNBC';
BC = 'SHBC';
switch BC
    case 'SHBC'
        noDomains = 1;
    case 'SSBC'
        noDomains = 2;
    case 'NNBC'
        noDomains = 3;
end
applyLoad = 'planeWave';
% applyLoad = 'pointCharge';
r_s = 6;

coreMethod = 'IGA';

varCol = setIhlenburgParameters(noDomains);
varCol{1}.meshFile = 'createNURBSmesh_EL';
varCol{1}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, max(2^(M-1)/8-1,0)];
if numel(varCol) > 1
    varCol{2}.refinement = @(M,t,t_fluid) [2^(M-1)-1, 2^(M-1)-1, max(round(t/t_fluid)*2^(M-1),0)];
end
if numel(varCol) > 2
    varCol{3}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, max(2^(M-2)-1,0)];
end
c_f = varCol{1}.c_f;   % Speed of sound in outer fluid
k = 2;
omega = k*c_f;
f = omega/(2*pi); 
parm = 1;


postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'method','formulation','M'};
postPlot(1).subFolderName 	= '../results/Ihlenburg3';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;

M = 5; % 5

N = 4;

alpha_s = pi;
beta_s = 0;

degree = 2;
calculateFarFieldPattern = 1;
calculateVolumeError = 1;
calculateSurfaceError = 1;
prePlot.plot2Dgeometry = 0;
prePlot.plot3Dgeometry = 0;
prePlot.plotControlPolygon = 0;

para.name                   = '';
para.plotResultsInParaview	= true;
para.plotTimeOscillation     = 0;
para.extraXiPts              = '3';
para.extraEtaPts             = '3'; 
para.extraZetaPts            = '3'; 

scaleForErrorPlot = 0;
plotMesh = 1;
loopParameters = {'M','method'};

collectIntoTasks