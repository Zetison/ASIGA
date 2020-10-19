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
runTasksInParallel = true;

% Upsilon = [22*sqrt(2)/3, 44*sqrt(3)/7]; % 3:1, 7:1
Upsilon = 22*sqrt(2)/3; % 3:1
Xi = [0,0,0,1,1,2,2,3,3,3]/3;

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
M = 4:5;
parm = 1;
alpha = 0;
beta = (-90:0.5:90)*pi/180;
prePlot.plot3Dgeometry = 0;
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plotArtificialBndry = false;        % Plot the artificial boundary for the IENSG method
computeCondNumber = 0;
calculateSurfaceError = 1;
calculateFarFieldPattern = false;     % Calculate far field pattern
refineThetaOnly = true;

% prePlot.resolution = [20,20,0];

degree = 5;

warning('off','NURBS:weights')
loopParameters = {'N','p_ie','s_ie','IElocSup', 'IEbasis','method', 'M', 'c_x'};

postPlot(1).xname       	= 's_ie';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType    	= 'semilogy';
postPlot(1).lineStyle   	= '*-';
postPlot(1).xLoopName     	= 's_ie';
postPlot(1).yScale          = 1/100;
% postPlot(1).legendEntries 	= {'method','formulation','M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_error(i_study);
% 
% postPlot(2) = postPlot(1);
% postPlot(2).yname        	= 'cond_number';
% postPlot(2).addCommands   	= [];

IElocSup = true;
s_ie = linspace(0.5,4,100);
maxN = 100;
% maxN = 40;
for p_ie = 1:5
    N = p_ie*2.^(1:floor(log(maxN/p_ie)/log(2)));
    N = 4*p_ie;
%     N = [3*p_ie,4*p_ie];
    % N = 3;
    collectIntoTasks
end

IEbasis	= 'Lagrange';
% IEbasis	= 'Chebyshev';
for p_ie = 3 %:5
    N = p_ie*2.^(1:floor(log(maxN/p_ie)/log(2)));
%     N = p_ie*round((2*p_ie:p_ie:100)/p_ie);
%     N = p_ie*2.^(1:7);
%     N = 2*p_ie:p_ie:20;
    N = 4*p_ie;
%     collectIntoTasks
end

IElocSup = false;
N = 1:19;
% N = 1:15;
% N = 4;
p_ie = NaN;
IEbasis	= 'Bernstein';
% collectIntoTasks
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
method = {'BA'};
useNeumanProj = 0;
solveForPtot = 0;
N = [N(1),N(end)];
formulation = {'SL2E'};
% collectIntoTasks

function addCommands_error(i_study)
if i_study == 6
    for p = 1:5
        error_safjan = importdata(['miscellaneous/refSolutions/Safjan2002tdi_Figure4_IElocSup1_p' num2str(p) '.csv']);
        loglog(error_safjan(:,1),error_safjan(:,2),'*','DisplayName',['Safjan, p = ' num2str(p)]);
    end
    error_safjan = importdata('miscellaneous/refSolutions/Safjan2002tdi_Figure4_IElocSup0.csv');
    loglog(error_safjan(:,1),error_safjan(:,2),'*','DisplayName','Safjan, Global multipole IE');
    legend('off');
    legend('show','Interpreter','latex');
end

