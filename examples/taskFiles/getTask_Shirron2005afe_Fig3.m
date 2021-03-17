function studies = getTask_Shirron2005afe_Fig3()
% This study is based on Hetmaniuk2012raa (Figure 7)
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271

counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
misc.scatteringCase = 'BI';
misc.model = 'SMS';  % Spherical shell
% misc.coreMethod = {'C0_IGA'};
misc.coreMethod = {'IGA'};
misc.applyLoad = 'planeWave';
% misc.applyLoad = 'radialPulsation';

% misc.method = {'IENSG'};
% misc.method = {'IE'};
misc.method = {'PML'};
BC = {'SHBC'};


calculateFarFieldPattern    = true;     % Calculate far field pattern
alpha_s = 0;                            % Aspect angle of incident wave
beta_s  = 0;                        % Elevation angle of incident wave
alpha   = linspace(0,pi,1000);                            % Aspect angles of observation points
beta = 0;   
r = 1;                            % radii for near-field evaluation.

a = 1;
N = 9; % 9
varCol{1} = struct('media', 'fluid', ...
                   't', 1, ...
                   'R1', a, ...
                   'R2', a, ...
                   'L', 4*a, ...
                   'c_f', 1500, ...
                   'rho', 1000);
varCol{1}.chimin = 2.9;
varCol{1}.chimax = 3.2;
varCol{1}.meshFile = 'createNURBSmesh_M3';
refineThetaOnly = true;
Xi = [0,0,0,1,1,2,2,3,3,3]/3;
% Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
switch misc.method{1}
    case 'PML'
        formulation = {'GSB'};
    %     formulation = {'STD'};
        varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-1)/8-1, 3*2^(M-1)/16-1];
    %     varCol{1}.refinement = @(M) [0, round(pi/0.5*1.5*10/2), 2, 6];
    %     varCol{1}.refinement = @(M) [0, 0, 4, 32];
    case {'IE','IENSG'}
        formulation = {'BGU'};
        varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-1)/8-1];
    %     varCol{1}.refinement = @(M) [0, 0, 4];
end

degree = 2;
M = 6; % 6
extraGP = [9-degree,0,0];    % extra quadrature points
% extraGP = [9-degree,9-degree,0];    % extra quadrature points
warning('off','NURBS:weights')

k = 10/a;
lambda = 2*pi/k;
f = k*varCol{1}.c_f/(2*pi);
omega = 2*pi*f;
% t_PML = 0.5*a;
r_a = 1.25*a;

pml.eps = 1e9*eps;      % choosing eps = eps yields machine precicion at Gamma_b, but requires more "radial" elements in the PML to resolve the rapid decay function
pml.sigmaType = 2;   	% sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n
pml.gamma = 5;          % parameter for sigmaType = 1
pml.t = 0.25*a;         % thickness of PML
pml.dirichlet = false;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)

IElocSup = 1;        % Toggle usage of radial shape functions in IE with local support
p_ie     = 4;          % Set polynomial order for radial shape functions
s_ie     = 2;          % Distrubution order for radial shape functions

parm = 1;
calculateFarFieldPattern = 1;
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;
% prePlot.colorFun = @(v) abs(norm2(v)-(r_a+t_PML));
% prePlot.resolution = [20,20,0];
prePlot.resolution = [20,0,0];
computeCondNumber = 0;

postPlot(1).xname           = 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'plot';
postPlot(1).lineStyle    	= '-';
postPlot(1).xScale       	= 180/pi;
postPlot(1).yScale       	= 1;
postPlot(1).legendEntries 	= {};
postPlot(1).subFolderName 	= '';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(i_study);

loopParameters = {'M','misc.method','BC'};
para.plotResultsInParaview	= 1;
para.extraXiPts              = '0';  % Extra visualization points in the xi-direction per element
para.extraEtaPts             = 'round(20/2^(M-1))';  % Extra visualization points in the eta-direction per element
para.extraZetaPts            = 'round(1/2^(M-1))';   % Extra visualization points in the zeta-direction per element
collectIntoTasks

misc.method = {'BEM'};
formulation = {'GBM'};
solveForPtot = true;
varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-1)/8-1];
% collectIntoTasks



function addCommands_(i_study)
if i_study == 1
    T = readtable('miscellaneous/refSolutions/SMS_M7_BEM_SHBC_TSVSalpha.txt', ...
                    'FileType','text', 'HeaderLines',1,'CommentStyle','%');
    x = T.Var1;
    y = T.Var2;
    plot(x,y,'DisplayName','Reference solution using GBM BEM')
    legend('off');
    legend('show','Interpreter','latex');
end


