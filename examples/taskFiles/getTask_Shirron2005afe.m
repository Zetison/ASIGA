function studies = getTask_Shirron2005afe()
% This study is based on Hetmaniuk2012raa (Figure 7)
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271

counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
hetmaniukCase = false; % evaluating solution at boundary not implemented

scatteringCase = 'BI';
model = 'S1';  % Spherical shell
coreMethod = {'IGA'};
applyLoad = 'planeWave';

% method = {'IENSG'};
% method = {'IE'};
method = {'PML'};
BC = {'SHBC'};

plotFarField = ~hetmaniukCase;
% plotFarField = true;     % If false, plots the near field instead

calculateFarFieldPattern    = true;     % Calculate far field pattern
alpha_s = 0;                            % Aspect angle of incident wave
beta_s  = -pi/2;                        % Elevation angle of incident wave
alpha   = 0;                            % Aspect angles of observation points
beta = linspace(-pi/2,pi/2,1000);   
r = 1;                            % radii for near-field evaluation.
extraGP = [7,0,0];    % extra quadrature points

a = 1; % Midsurface radius
varCol{1} = struct('media', 'fluid', ...
                  'R_i', a, ...
                  't',   a, ...
                  'c_f', 1524, ...
                  'rho', 1000);
N = 4; % 9
varCol{1}.meshFile = 'createNURBSmesh_EL';
Xi = [0,0,0,1,1,2,2,3,3,3]/3;
if strcmp(method{1},'PML')
%     formulation = {'GSB'};
    formulation = {'STD'};
    varCol{1}.refinement = @(M) [0, 2^(M-1)-1, max(2^(M-1)/8-1,10)];
else
    formulation = {'BGU'};
    varCol{1}.refinement = @(M) [0, 2^(M-1)-1, max(2^(M-1)/8-1,10-N+1)];
end

degree = 2;
M = 7; % 6

k = 1/a;
k = 10/a;
lambda = 2*pi/k;
f = k*varCol{1}.c_f/(2*pi);
omega = 2*pi*f;
r_PML = 1.25*a;
r_PML = a;
r_a = 1.5*a;
gamma = 5;
% gamma = 10;

parm = 1;
calculateSurfaceError = 1;
calculateVolumeError  = 0;
calculateFarFieldPattern = 1;
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;
% prePlot.colorFun = @(v) abs(norm2(v)-1);
prePlot.resolution = [20,20,0];
computeCondNumber = 0;

postPlot(1).xname           = 'M';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'semilogy';
postPlot(1).lineStyle    	= '-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).legendEntries 	= {};
postPlot(1).subFolderName 	= '';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).xLoopName     	= 'M';

postPlot(2) = postPlot(1);
postPlot(2).xname           = 'beta';
postPlot(2).yname           = 'TS';
postPlot(2).axisType      	= 'plot';
postPlot(2).noXLoopPrms   	= 0;

postPlot(3) = postPlot(2);
postPlot(3).axisType = 'semilogy';
postPlot(3).yname = 'error_p';

loopParameters = {'M','method','BC'};
para.plotResultsInParaview	= 0;
para.extraXiPts              = '30';  % Extra visualization points in the xi-direction per element
para.extraEtaPts             = 'round(20/2^(M-1))';  % Extra visualization points in the eta-direction per element
para.extraZetaPts            = 'round(1/2^(M-1))';   % Extra visualization points in the zeta-direction per element
collectIntoTasks
