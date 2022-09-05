function studies = getTask_Fillinger()
% This study is based on Figure 14 page 292 in Venas2019asi
% Venas2019asi is available at http://hdl.handle.net/11250/2640443

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

saveStudies = false;

misc.scatteringCase = 'Sweep';

misc.model = 'S1';  % Spherical shell

% misc.coreMethod = {'IGA', 'XI'};
misc.coreMethod = {'IGA'};
misc.method = {'BEM'};
if strcmp(misc.method, 'BEM')
    misc.formulation = {'CCBIE', 'CBM', 'CHBIE'};
end

varCol = setS1Parameters('double',1);
msh.meshFile = 'createNURBSmesh_EL';

% Eigenfrequencies from Zheng2015itb (note that they are multiplied by two, to
% compensate for the half size of the sphere)
eigenValues = [3.141592653589794 % pi
                   6.283185307179587 % 2*pi
                   9.424777960769379 % 3*pi
                   4.493409457909065
                   7.725251836937707
                   5.763459196894550
                   9.095011330476353
                   6.987932000500519
                   8.182561452571243
                   9.355812111042747]';

prePlot.plot3Dgeometry = 0;
prePlot.abortAfterPlotting = 1;       % Abort simulation after pre plotting
k = 10.^linspace(-1,4,1000); % 10.^linspace(-1,2,1000)
% k = 10.^linspace(1,2,20); % 10.^linspace(-1,2,1000)
% k = sort([eigenValues, k]);
% eigenValues = eigenValues*1500/(2*pi);
c_f = 1500;
misc.omega = k*c_f;

ffp.alpha_s = 0;
ffp.beta_s = pi/2;   
ffp.alpha = ffp.alpha_s;
ffp.beta = ffp.beta_s;  

msh.M = 1:3;
err.calculateSurfaceError = 1;

postPlot(1).xname        	= 'k';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= 1;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'semilogx';
postPlot(1).lineStyle    	= '-';

% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
misc.method = {'MFS'};
misc.coreMethod = 'IGA';
misc.formulation = 'PS';
msh.M = 5;
msh.degree = 2;
err.calculateSurfaceError = 0;
loopParameters = {'misc.omega','msh.M','misc.method'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KDT simulation
misc.method = {'KDT'};
misc.formulation = 'MS1';
msh.degree = 4;
err.calculateSurfaceError = 0;
loopParameters = {'misc.method','msh.M','msh.parm','misc.coreMethod'};
misc.coreMethod = {'linear_FEM'};
% misc.coreMethod = {'IGA'};
msh.M = [6,8,10];
msh.M = 6;
msh.parm = 2;
collectIntoTasks


misc.coreMethod = {'IGA'};
msh.M = 6;
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT simulation
misc.method = {'RT'};
misc.formulation = '';
msh.M = 3;
msh.degree = 2;
err.calculateSurfaceError = 0;
misc.plotFarField = 1;
misc.applyLoad = 'planeWave';
rt.N = 3:6;
msh.parm = 1;

loopParameters = {'misc.method','msh.M','msh.parm','rt.N','misc.coreMethod'};
% collectIntoTasks
