%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 18 in Venas2018iao
% Venas2018iao is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'IL';

method = {'IE'};
formulation = 'BGU';

BC = 'NNBC';

coreMethod = 'IGA';

varCol = setIhlenburgParameters(3);
varCol{1}.meshFile = 'createNURBSmesh_EL';
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

M = 5;

N = 6;

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
para.extraXiPts              = '3';
para.extraEtaPts             = '3'; 
para.extraZetaPts            = '3'; 

scaleForErrorPlot = 0;
plotMesh = 1;
loopParameters = {'M','method'};

collectIntoTasks