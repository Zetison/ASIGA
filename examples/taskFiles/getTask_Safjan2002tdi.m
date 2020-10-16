function studies = getTask_Safjan2002tdi()
% This study is based on Simpson2014aib and correspond to Figure 6 in Venas2018iao
% Simpson2014aib is available at https://doi.org/10.1016/j.cma.2013.10.026
% Venas2018iao   is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

%% IE simulation
scatteringCase = 'BI';
model = 'Safjan2002tdi'; % Simpson sphere
method = {'IENSG'};
% method = {'IE'};
IEbasis	= 'Chebyshev';
% IEbasis	= 'Lagrange';
% IEbasis	= {'Chebyshev','Bernstein','Lagrange'};
BC = 'NBC';
formulation = {'PGU','PGC','BGU','BGC'};
formulation = {'PGC','BGC'};
formulation = {'PGC'};
coreMethod = 'IGA';
runTasksInParallel = 0;
progressBars = false;        % Show progress bars for building system matrices
applyLoad = 'Safjan';

% Upsilon = [22*sqrt(2)/3, 44*sqrt(3)/7]; % 3:1, 7:1
Upsilon = 22*sqrt(2)/3; % 3:1

c_z = 11;
c_x = sqrt(c_z^2-Upsilon.^2);

layer{1} = struct('media', 'fluid', ...
                  'c_f', 1500, ...
                  'rho', 1000);
              
varCol = setSSParameters(1);
varCol{1}.meshFile = 'createNURBSmesh_EL';
c_f = varCol{1}.c_f;   % Speed of sound in outer fluid
k = 1;                 % Wave number for Simpson2014aib
omega = c_f*k;         % Angular frequency
f = omega/(2*pi);      % Frequency

M = 1:7;
M = 5;
parm = 1;
alpha = 0;
beta = (-90:0.5:90)*pi/180;
prePlot.plot3Dgeometry = 0;
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plotArtificialBndry = false;        % Plot the artificial boundary for the IENSG method
computeCondNumber = 1;
calculateSurfaceError = 1;
calculateFarFieldPattern = false;     % Calculate far field pattern
% prePlot.resolution = [20,20,0];

degree = 5;

warning('off','NURBS:weights')
loopParameters = {'N','p_ie','IElocSup', 'c_x'};

postPlot(1).xname       	= 'N';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType    	= 'loglog';
postPlot(1).lineStyle   	= '*-';
postPlot(1).xLoopName     	= 'N';
postPlot(1).yScale          = 1/100;
% postPlot(1).legendEntries 	= {'method','formulation','M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).addCommands   	= [];

postPlot(2) = postPlot(1);
postPlot(2).yname        	= 'cond_number';

% postPlot(3) = postPlot(1);
% postPlot(3).noXLoopPrms     = 0;
% postPlot(3).lineStyle       = '-';
% postPlot(3).xname           = 'beta';
% postPlot(3).yname           = 'error_p';
% postPlot(3).axisType        = 'semilogy';
% postPlot(3).xScale          = 180/pi;

% postPlot(4) = postPlot(3);
% postPlot(4).yname           = 'abs_p';
% postPlot(4).axisType        = 'plot';

IElocSup = true;
for p_ie = 3:5
    N = 2*p_ie:11;
    % N = 3;
%     collectIntoTasks
end

IEbasis	= 'Lagrange';
for p_ie = 3:5
    N = 2*p_ie:p_ie:11;
    N = 4*p_ie;
%     collectIntoTasks
end

IElocSup = false;
N = 1:11;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
method = {'BA'};
useNeumanProj = 0;
solveForPtot = 0;
N = [N(1),N(end)];
formulation = {'SL2E'};
% collectIntoTasks

