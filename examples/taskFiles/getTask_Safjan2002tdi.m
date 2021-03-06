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
coreMethod = 'IGA';
runTasksInParallel = 0;
progressBars = false;        % Show progress bars for building system matrices
applyLoad = 'Safjan';

Upsilon = [22*sqrt(2)/3, 44*sqrt(3)/7]; % 3:1, 7:1
Upsilon = 22*sqrt(2)/3; % 3:1
% Upsilon = 44*sqrt(3)/7; % 7:1
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
M = 6; %8
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

degree = 5;

warning('off','NURBS:weights')
loopParameters = {'N','p_ie','s_ie','IElocSup', 'IEbasis','method', 'formulation', 'c_x'};

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
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_error();
% 
% postPlot(2) = postPlot(1);
% postPlot(2).yname        	= 'cond_number';
% postPlot(2).addCommands   	= [];

runTasksInParallel = 0;
IElocSup = true;
s_ie = [1,2];
s_ie = 1;
N_arr = {[2,10,20,40,90],...
         [4,20,30,40,80,100],...
         [6,30,40,60,90],...
         [8,20,40,60,80],...
         [10,25,35,50,60,75,100]};
% formulation = {'PGC','BGU'};
% formulation = {'BGU'};
formulation = {'PGU'};
% formulation = {'WBGC','PGC','WBGU','PGU'};
maxN = 100;
% maxN = 40;
for p_ie = 3 %1:5
%     N = p_ie*2.^(1:floor(log(maxN/p_ie)/log(2)));
    N = N_arr{p_ie};
%     N = 4*p_ie;
%     N = [3*p_ie,4*p_ie];
%     N = 2*p_ie:p_ie:4*p_ie;
    % N = 3;
    collectIntoTasks
end

IEbasis	= 'Lagrange';

for p_ie = 3 %1:5
%     N = p_ie*2.^(1:floor(log(maxN/p_ie)/log(2)));
    N = N_arr{p_ie};
%     N = p_ie*round((2*p_ie:p_ie:100)/p_ie);
%     N = p_ie*2.^(1:7);
%     N = 2*p_ie:p_ie:20;
%     N = 4*p_ie;
%     collectIntoTasks
end

IElocSup = false;
s_ie = NaN;
N = 1:19;
N = 1:9;
% N = 4;
p_ie = NaN;
IEbasis	= 'Chebyshev';
IEbasis	= 'Lagrange';
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
method = {'BA'};
useNeumanProj = 0;
solveForPtot = 0;
N = [1,maxN];
formulation = {'SL2E'};
postPlot(1).lineStyle = '--';
collectIntoTasks

function addCommands_error()
for ar = [3,7]
    for formulation = {'PGC','BGU'}
        for s_ie = 1:2
            for p_ie = 1:5
                try
                    error_safjan = importdata(['miscellaneous/refSolutions/Safjan2002tdi_IElocSup1_' formulation{1} '_s_ie' num2str(s_ie) '_p_ie' num2str(p_ie) '_ar' num2str(ar) ':1_surfaceError.csv']);
                    loglog(error_safjan(:,1),error_safjan(:,2),'-','DisplayName',['Safjan, s = ' num2str(s_ie) ', p = ' num2str(p_ie) ', ' formulation{1} ' IE, ' num2str(ar) ':1']);
                end
            end
        end
        error_safjan = importdata(['miscellaneous/refSolutions/Safjan2002tdi_IElocSup0_' formulation{1} '_s_ieNaN_p_ieNaN_ar' num2str(ar) ':1_surfaceError.csv']);
        loglog(error_safjan(:,1),error_safjan(:,2),'-','DisplayName',['Safjan, Global multipole ' formulation{1} ' IE, ' num2str(ar) ':1']);
    end
end
legend('off');
legend('show','Interpreter','latex');

