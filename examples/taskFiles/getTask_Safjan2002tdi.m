function studies = getTask_Safjan2002tdi()
% This study is based on Simpson2014aib and correspond to Figure 6 in Venas2018iao
% Simpson2014aib is available at https://doi.org/10.1016/j.cma.2013.10.026
% Venas2018iao   is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

%% IE simulation
misc.scatteringCase = 'BI';
misc.model = 'Safjan2002tdi'; % Simpson sphere
misc.method = {'IENSG'};
% misc.method = {'IE'};
% iem.IEbasis	= 'Chebyshev';
iem.IEbasis	= 'Lagrange';
% iem.IEbasis	= {'Chebyshev','Bernstein','Lagrange'};
misc.coreMethod = 'IGA';
runTasksInParallel = 0;
progressBars = false;        % Show progress bars for building system matrices
misc.applyLoad = 'Safjan';
% misc.applyLoad = 'Safjan10';
% misc.applyLoad = 'planeWave';
switch misc.applyLoad
    case {'Safjan','Safjan10'}
        misc.BC = 'NBC';
    case 'planeWave'
        misc.BC = 'SHBC';
end
misc.checkNURBSweightsCompatibility = false;

Upsilon = [22*sqrt(2)/3, 44*sqrt(3)/7]; % 3:1, 7:1
Upsilon = 22*sqrt(2)/3; % 3:1
% Upsilon = 44*sqrt(3)/7; % 7:1
msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;

c_z = 11;
c_x = sqrt(c_z^2-Upsilon.^2);
% c_x = c_z;
              
varCol = setSSParameters(1);
c_f = varCol{1}.c_f;   % Speed of sound in outer fluid
k = 1;                 % Wave number for Simpson2014aib
misc.omega = c_f*k;    % Angular frequency

varCol{1} = struct('media', 'fluid', ...
                   'c_x',   c_x, ...
                   'c_z',   c_z, ...
                   'c_f', 1500, ...
                   'rho', 1000);
msh.meshFile = 'createNURBSmesh_EL';
c_f = varCol{1}.c_f;   % Speed of sound in outer fluid
k = 1;                 % Wave number for Simpson2014aib
misc.omega = c_f*k;         % Angular frequency

msh.M = 1:7;
msh.M = 5; %8
msh.parm = 1;
msh.explodeNURBS = 0;
msh.autoRefine = 0;   % Use automatic refining algorithm based on h_max = refLength/2^(M-1) where M is the mesh number and refLength is the largest element edge in the coarse mesh
ffp.alpha = 0;
ffp.beta = (-90:0.5:90)*pi/180;
ffp.alpha_s = 0;
ffp.beta_s = -pi/2;
prePlot.plot3Dgeometry = 0;
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plotArtificialBndry = false;        % Plot the artificial boundary for the IENSG misc.method
misc.computeCondNumber = 1;
err.calculateSurfaceError = 1;
ffp.calculateFarFieldPattern = 0;     % Calculate far field pattern
varCol{1}.refinement = @(M) [0, 2^(M-1)-1, max(2^(M-1)/8-1,0)];
msh.refineThetaOnly = true;
ffp.extraGP = [50,0,0];
misc.extraGP = [10,0,0];

msh.degree = 4;

warning('off','NURBS:weights')
loopParameters = {'iem.N','iem.p_ie','iem.s_ie','iem.IElocSup', 'iem.IEbasis','misc.method', 'misc.formulation', 'varCol{1}.c_x'};

postPlot(1).xname       	= 'iem.N';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType    	= 'loglog';
postPlot(1).lineStyle   	= '*-';
postPlot(1).xLoopName     	= 'iem.N';
postPlot(1).yScale          = 1/100;
% postPlot(1).legendEntries 	= {'misc.method','formulation','M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_error();

if misc.computeCondNumber
    postPlot(2) = postPlot(1);
    postPlot(2).yname        	= 'cond_number';
    postPlot(2).addCommands   	= [];
end
if ffp.calculateFarFieldPattern
    postPlot(2) = postPlot(1);
    postPlot(2).yname        	= 'TS';
    postPlot(2).xname       	= 'ffp.beta';
    postPlot(2).noXLoopPrms   	= 0;
    postPlot(2).addCommands   	= [];
    postPlot(2).xScale          = 180/pi;
    postPlot(2).yScale          = 1;
    postPlot(2).axisType    	= 'plot';
    postPlot(2).lineStyle   	= '-';
end

runTasksInParallel = 0;
iem.IElocSup = 1;
iem.s_ie = [1,2];
iem.s_ie = 1;
N_arr = {[2,10,20,40,90],...
         [4,20,30,40,80,100],...
         [6,30,40,60,90],...
         [8,20,40,60,80],...
         [10,25,35,50,60,75,100]};
% misc.formulation = {'PGC','BGU'};
% misc.formulation = {'BGU'};
misc.formulation = {'PGU'};
% misc.formulation = {'WBGC','PGC','WBGU','PGU'};
maxN = 100;
% maxN = 40;
for p_ie = 3 %1:5
    iem.p_ie = p_ie;
%     iem.N = p_ie*2.^(1:floor(log(maxN/p_ie)/log(2)));
    iem.N = N_arr{p_ie};
%     iem.N = 4*p_ie;
%     iem.N = [3*p_ie,4*p_ie];
%     iem.N = 2*p_ie:p_ie:4*p_ie;
    % iem.N = 3;
%     collectIntoTasks
end

iem.IEbasis	= 'Lagrange';

for p_ie = 3 %1:5
%     iem.N = p_ie*2.^(1:floor(log(maxN/p_ie)/log(2)));
    iem.N = N_arr{p_ie};
%     iem.N = p_ie*round((2*p_ie:p_ie:100)/p_ie);
%     iem.N = p_ie*2.^(1:7);
%     iem.N = 2*p_ie:p_ie:20;
%     iem.N = 4*p_ie;
%     collectIntoTasks
end

iem.IElocSup = false;
iem.s_ie = NaN;
% iem.N = 1:19;
iem.N = 1:12;
% iem.N = 4;
iem.p_ie = NaN;
% iem.IEbasis	= 'Chebyshev';
misc.formulation = {'PGC'};
% misc.formulation = {'BGU'};
iem.IEbasis	= 'Lagrange';
% iem.IEbasis	= 'Chebyshev';
iem.boundaryMethod = true;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
misc.method = {'BA'};
misc.useNeumanProj = 0;
misc.solveForPtot = 0;
iem.N = [1,maxN];
misc.formulation = {'SL2E'};
postPlot(1).lineStyle = '--';
collectIntoTasks

function addCommands_error()
for ar = 3 %[3,7]
    for formulation = {'PGC','BGU'}
        for s_ie = 1:2
            for p_ie = 1:5
                try
                    error_safjan = importdata(['miscellaneous/refSolutions/Safjan2002tdi_IElocSup1_' formulation{1} '_s_ie' num2str(s_ie) '_p_ie' num2str(p_ie) '_ar' num2str(ar) '_1_surfaceError.csv']);
                    loglog(error_safjan(:,1),error_safjan(:,2),'-','DisplayName',['Safjan, s = ' num2str(s_ie) ', p = ' num2str(p_ie) ', ' formulation{1} ' IE, ' num2str(ar) '\_1']);
                end
            end
        end
        error_safjan = importdata(['miscellaneous/refSolutions/Safjan2002tdi_IElocSup0_' formulation{1} '_s_ieNaN_p_ieNaN_ar' num2str(ar) '_1_surfaceError.csv']);
        loglog(error_safjan(:,1),error_safjan(:,2),'-','DisplayName',['Safjan, Global multipole ' formulation{1} ' IE, ' num2str(ar) '\_1']);
    end
end
legend('off');
legend('show','Interpreter','latex');
% legend('show');

