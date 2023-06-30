function studies = getTask_CSLP()
% This study is based on Hetmaniuk2012raa (Figure 7)
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271

counter = 1;
studies = cell(0,1);
getDefaultTaskValues
saveStudies        = false;       % save ASIGA-struct into a .mat file
appyCommandsAt = 1;
noCoresToUse       = 4;

%% IE simulation
% BC = 'NNBC';
misc.model = 'S1';  % Spherical shell
misc.coreMethod = {'IGA'};
misc.applyLoad = 'planeWave';

sol.preconditioner = {'CSLP','diag'};
sol.solver = {'gmres','LU'};
connectedParameters = {{'sol.solver','sol.preconditioner'},{'msh.M','misc.omega','pml.gamma'}};

BCs = {'SHBC'};
% BCs = {'SSBC'};
% BCs = {'NNBC'};
% BCs = {'SHBC','SSBC','NNBC'};

% sol.solver          = 'gmres';  % 'LU', 'gmres', 'cgs', 'bicgstab', 'bicgstabl', 'lsqr', 'bicg'
% sol.preconditioner  = 'ilu';	% 'ilu', 'SSOR', 'diag'

warning('off','NURBS:weights')
ffp.alpha_s = 240*pi/180;                            % Aspect angle of incident wave
ffp.beta_s  = 30*pi/180;                        % Elevation angle of incident wave

msh.parm = 1;
msh.explodeNURBS = 0;   % Create patches from all C^0 interfaces
msh.refineThetaOnly = 0; % only for msh.parm = 1
err.calculateSurfaceError = 1;
err.calculateVolumeError  = 0;
misc.calculateFarFieldPattern = 1;
misc.checkNURBSweightsCompatibility = false;
misc.preProcessOnly = 0;

prePlot.plot3Dgeometry = 0;
prePlot.view                = [23,18];     % Set view angle [azimuth,elevation]
prePlot.plotGeometryInfo    = false;      % Plot domain boundaries (i.e. Gamma, Gamma_a, Neumann, Dirichlet, ...)
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plotFullDomain      = 1;
prePlot.plotControlPolygon  = 0;
% prePlot.plotSubsets         = {'xz'};
prePlot.plotSubsets         = {}; 
prePlot.view                = [0,0];     % Set view angle [azimuth,elevation]
% prePlot.colorFun = @(v) abs(norm2(v)-1);
prePlot.resolution = [100,40,0];
% prePlot.resolution = [400,200,0];
% prePlot.resolution = [100,0,0];
if prePlot.plotFullDomain
    prePlot.format = '-png';      % Use '-png' or '-pdf' (for vector graphics)
else
    if strcmp(prePlot.plotSubsets{1},'xz')
        prePlot.format = '-pdf';      % Use '-png' or '-pdf' (for vector graphics)
    else
        prePlot.format = '-png';      % Use '-png' or '-pdf' (for vector graphics)
    end
end

for i = 1:numel(BCs)
    misc.BC = BCs{i};
%     misc.method = {'IENSG'};
%     misc.method = {'IE'};
    misc.method = {'PML'};

    postPlot(1).xname           = 'sol.beta_CSLP';
%     postPlot(1).xname           = 'surfDofs';
    postPlot(1).yname        	= 'surfaceError';
    postPlot(1).plotResults  	= true;
    postPlot(1).printResults 	= true;
%     postPlot(1).axisType      	= 'semilogy';
%     postPlot(1).axisType      	= 'loglog';
    postPlot(1).axisType      	= 'semilogy';
    postPlot(1).lineStyle    	= '*-';
    postPlot(1).xScale       	= 1;
    postPlot(1).yScale       	= 1;
    postPlot(1).legendEntries 	= {};
    postPlot(1).subFolderName 	= '';
    postPlot(1).fileDataHeaderX	= [];
    postPlot(1).noXLoopPrms   	= 0;
    postPlot(1).xLoopName       = NaN;
    postPlot(1).plotResults     = true;
    postPlot(1).printResults    = true;
    postPlot(1).noXLoopPrms   	= 1;
%     postPlot(1).xLoopName     	= 'msh.M';
    postPlot(1).xLoopName     	= 'sol.beta_CSLP';
    postPlot(2)                 = postPlot(1);
    postPlot(2).yname        	= 'timeSolveSystem';
%     postPlot(2).axisType      	= 'loglog';
%     postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(study);
    switch misc.BC
        case {'SHBC','NBC'}
            noDomains = 1;
        case 'SSBC'
            noDomains = 2;
        case 'NNBC'
            noDomains = 3;
    end
    varCol = setHetmaniukParameters(noDomains);
    msh.meshFile = 'createNURBSmesh_EL';
    if msh.refineThetaOnly
        msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
    else
        msh.Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
    end
    msh.degree = 2;
    manuelRefinement = 1;
    misc.formulation = {'GSB'};
    if msh.refineThetaOnly
        misc.extraGP = [9-msh.degree(1),0,0];    % extra quadrature points
        ffp.extraGP = [50,0,0];    % extra quadrature points
    end

    msh.M = 5:7; % 5:7
%     msh.M = 1:2;
    misc.symmetric = 0;

    iem.N = 3; % 9
    iem.p_ie = 5;
    iem.s_ie = 2;
    iem.IElocSup = 1;        % Toggle usage of radial shape functions in IE with local support
    
    pml.eps = 1e9*eps;      % choosing eps = eps yields machine precicion at Gamma_b, but requires more "radial" elements in the PML to resolve the rapid decay function
    pml.sigmaType = 3;  % sigmaType = 1: sigma(xi) = xi*(exp(gamma*xi)-1), sigmaType = 2: sigma(xi) = C*xi^n
    pml.dirichlet = 1;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
    
    k = 5*2.^(msh.M-5);
    c_f = varCol{1}.c_f;
    f = k*c_f/(2*pi);
    omega = 2*pi*f;
    misc.omega = omega;

    pml.t = 0.2*varCol{1}.R;         % thickness of PML
    pml.gamma = 1./(k*pml.t);
    misc.r_a = 1.2*varCol{1}.R;

    misc.storeFullVarCol = false;
    sol.droptol         = [1e-2,1e-4];     % parameter for incomplete lu factorization with threshold and pivoting (ilutp)
    sol.beta_CSLP       = [linspace(0,0.1,11),0.5]; 
    
    misc.scatteringCase = 'BI';
    loopParameters = {'msh.M','misc.method','misc.coreMethod','misc.BC','misc.omega','sol.solver','sol.droptol','sol.beta_CSLP'};
    collectIntoTasks

    %% Run BA sweep
    misc.omega = omega;
%     misc.omega = misc.omega(1);
    para.plotResultsInParaview = 0;
    misc.method = {'BA'};
    misc.formulation = {'SL2E'};
%     misc.formulation = {'VL2E'};
%     collectIntoTasks
end

% function addCommands_(study)
% for i_task = 1:numel(study.tasks)
%     task = study.tasks(i_task).task;
%     c_f = task.varCol{1}.c_f;
% 
%     if isfield(task.rom,'history')
%         history = task.rom.history;
%         for i = 1:numel(history)
%             figure(20+i)
%             omega = task.rom.history(i).omega;
%             residual = task.rom.history(i).residual;
%             hold on
%         
%             omega_P = history(i).omega_P;
%             J_P = history(i).J_P;
%             semilogy(omega_P/c_f,1e-15*ones(size(omega_P)),'o','color','magenta','DisplayName','Interpolation points')
%             for j = 1:numel(omega_P)
%                 text((omega_P(j)+omega_P(end)/200)/c_f,1e-15,num2str(J_P(j)),'color','magenta')
%             end
%     
%             semilogy([omega_P(1),omega_P(end)]/c_f,task.rom.tolerance*ones(1,2),'red','DisplayName','Tolerance')
%     
%             semilogy(omega/c_f,residual,'*','color','black','DisplayName','Residual check points')
%     
%             omega_T_new = history(i).omega_T_new;
%             [~,idx] = ismember(omega,omega_T_new);
%             semilogy(omega_T_new/c_f,residual(logical(idx)),'o','color','cyan','DisplayName','New interpolation point')
%     
%             ylim([5e-16,10])
%             set(gca,'yscale','log')
%             legend show
%             savefig([task.resultsFolder, '/', task.saveName '_adaptiveROM_iter' num2str(i)])
%         end
%         legend off
%         if isfield(task.rom.history(end),'residualFine')
%             residual = task.rom.history(end).residualFine;
%             omegaFine = task.rom.history(end).omegaFine;
%             semilogy(omegaFine/c_f,residual,'blue','DisplayName','Final relative residual')
%         end
%         if isfield(task.rom.history(end),'relError')
%             relError = task.rom.history(end).relError;
%             semilogy(omegaFine/c_f,relError,'green','DisplayName','Final relative error')
%         end
%         set(gca,'yscale','log')
%         ylim([5e-16,10])
%         ylabel('Relative error/residual')
%         xlabel('Wavenumber')
%         legend show
%         savefig([task.resultsFolder, '/', task.saveName '_adaptiveROM_iter' num2str(i)])
%     end
% end
% 



