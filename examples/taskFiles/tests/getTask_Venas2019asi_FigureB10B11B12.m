function studies = getTask_Venas2019asi_FigureB10B11B12(M_0)
% This study is based on Figure B.10, B.11 and B.12 in Venas2019asi
% Venas2019asi is available at http://hdl.handle.net/11250/2640443

if nargin < 1
    M_0 = 6; 
end

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'BI';

misc.model = 'S1';  % Spherical shell

misc.coreMethod = 'IGA';
varCol = setS1Parameters('double',1);
msh.meshFile = 'createNURBSmesh_EL';


f = 20e3;
misc.omega = 2*pi*f;

ffp.alpha_s = 0;
ffp.beta_s = 0;   

postPlot(1).xname        	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= 0;
postPlot(1).printResults 	= 0;
postPlot(1).axisType      	= 'plot';
postPlot(1).lineStyle    	= '-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).xLoopName     	= 'rt.N';

postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;

postPlot(2) = postPlot(1);
postPlot(2).noXLoopPrms = 0;
postPlot(2).lineStyle = '-';
postPlot(2).xname = 'alpha';
postPlot(2).yname = 'error_pAbs';
postPlot(2).axisType = 'semilogy';
postPlot(2).xScale = 180/pi;

postPlot(3) = postPlot(2);
postPlot(3).yname = 'error_p';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT simulation
misc.method = {'RT'};
misc.formulation = '';
msh.M = 3;
prePlot.plot3Dgeometry = 0;
msh.degree = 2;
err.calculateSurfaceError = 0;
misc.computeCondNumber = false;
ffp.plotFarField = 1;
misc.applyLoad = 'planeWave';
msh.parm = 1;
rt.N = M_0+1;
warning('off','RT:limitations')

loopParameters = {'rt.N'};
collectIntoTasks
