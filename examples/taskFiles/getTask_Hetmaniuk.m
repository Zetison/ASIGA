function studies = getTask_Hetmaniuk()
% This study is based on Hetmaniuk2012raa (Figure 7)
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271

counter = 1;
studies = cell(0,1);
getDefaultTaskValues
saveStudies        = false;       % save ASIGA-struct into a .mat file
appyCommandsAt = 1;

%% IE simulation
hetmaniukCase = 1;

% BC = 'NNBC';
misc.model = 'Hetmaniuk2012raa';  % Spherical shell
misc.coreMethod = {'IGA','C0_IGA'};
% misc.coreMethod = {'C0_IGA'};
% misc.coreMethod = {'IGA'};
misc.coreMethod = {'hp_FEM'};
% misc.applyLoad = 'pointPulsation';
% misc.applyLoad = 'pointCharge';
misc.applyLoad = 'planeWave';

BCs = {'SHBC'};
% BCs = {'SSBC'};
% BCs = {'NNBC'};
% BCs = {'SHBC','SSBC','NNBC'};
if strcmp(misc.applyLoad,'pointPulsation')
    BCs = {'NBC'};
end
ffp.plotFarField = ~hetmaniukCase;
% plotFarField = true;     % If false, plots the near field instead

ffp.calculateFarFieldPattern    = true;     % Calculate far field pattern
ffp.alpha_s = 0;                            % Aspect angle of incident wave
ffp.beta_s  = -pi/2;                        % Elevation angle of incident wave
ffp.alpha   = 0;                            % Aspect angles of observation points
ffp.beta = pi/2;   
ffp.r = 1;                            % radii for near-field evaluation.
ffp.splineBasedNFPcalc = true;
warning('off','NURBS:weights')

misc.r_s = 3;

msh.parm = 2;
msh.explodeNURBS = 0;   % Create patches from all C^0 interfaces
msh.refineThetaOnly = 0; % only for msh.parm = 1
err.calculateSurfaceError = 1;
err.calculateVolumeError  = 0;
misc.calculateFarFieldPattern = 1;
misc.checkNURBSweightsCompatibility = false;

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

misc.computeCondNumber = 0;

for i = 1:numel(BCs)
    misc.BC = BCs{i};
%     misc.method = {'IENSG'};
%     misc.method = {'IE'};
    misc.method = {'PML'};

    postPlot(1).xname           = 'varCol{1}.k';
    postPlot(1).yname        	= 'surfaceError';
    postPlot(1).plotResults  	= true;
    postPlot(1).printResults 	= true;
    postPlot(1).axisType      	= 'semilogy';
    postPlot(1).lineStyle    	= '-';
    postPlot(1).xScale       	= 1;
    postPlot(1).yScale       	= 1;
    postPlot(1).legendEntries 	= {};
    postPlot(1).subFolderName 	= '';
    postPlot(1).fileDataHeaderX	= [];
    postPlot(1).noXLoopPrms   	= 0;
    postPlot(1).xLoopName       = NaN;
    postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(study);

    postPlot(2) = postPlot(1);
    if hetmaniukCase
        postPlot(2).yname = 'p_Re';
    else
        postPlot(2).yname = 'TS';
    end
    postPlot(2).axisType      	= 'plot';

    postPlot(3) = postPlot(2);
    postPlot(3).axisType = 'semilogy';
    postPlot(3).yname = 'error_p';

    postPlot(4) = postPlot(1);
    postPlot(5) = postPlot(2);
    postPlot(6) = postPlot(3);
    for j = 4:6
        postPlot(j).xname = 'f';
    end
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
    msh.degree = 3;
    manuelRefinement = 1;
    if strcmp(misc.method{1},'PML')
        misc.formulation = {'GSB'};
    else
        misc.formulation = {'BGC'};
    end
    if manuelRefinement
        if msh.parm == 1
            if msh.refineThetaOnly
                varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-4)-1];
            else
                varCol{1}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, 2^(M-4)-1];
            end
        else
            varCol{1}.refinement = @(M) [3*2^(M-3)-1, 3*2^(M-3)-1, 2^(M-1)/8-1];
        end
    end
    if manuelRefinement
        pml.refinement = @(M) max(2^(M-4)-1,3);
    end
    misc.extraGP = [9-msh.degree(1),0,0];    % extra quadrature points
    ffp.extraGP = [50,0,0];    % extra quadrature points
    if noDomains > 1
        if manuelRefinement
            if msh.parm == 1
                if msh.refineThetaOnly
                    varCol{2}.refinement = @(M,t,t_fluid) [0, 2^(M-1)-1, max(round(t/t_fluid)*2^(M-1),0)];
    %                 varCol{2}.refinement = @(M,t,t_fluid) [0, 2^(M-1)-1, max(round(t/t_fluid*2^(M-3)),0)];
                else
                    varCol{2}.refinement = @(M,t,t_fluid) [2^(M-1)-1, 2^(M-1)-1, max(round(t/t_fluid*2^(M-4)),0)];
                end
            else
                varCol{2}.refinement = @(M,t,t_fluid) [3*2^(M-3)-1, 3*2^(M-3)-1, max(round(t/t_fluid*2^(M-4)),0)];
            end
        end
    end

    msh.M = 6; % 7 (6 if hetmaniukCase and manuelRefinement)
    rom.basisROM = {'Pade','Taylor','DGP','Hermite','Bernstein'};  % do not put basisROMcell in loopParameters (this is done automatically)
    rom.basisROM = {'Pade','DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
    rom.basisROM = {'DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
    rom.adaptiveROM = 1;
    rom.computeROMresidualFine = 1;
    rom.computeROMerror = 1;
    rom.J_max = 20;
%     sol.preconditioner = 'none';
    misc.symmetric = 0;

    iem.N = 3; % 9
    iem.p_ie = 5;
    iem.s_ie = 2;
    iem.IElocSup = 1;        % Toggle usage of radial shape functions in IE with local support
    
    pml.eps = 1e9*eps;      % choosing eps = eps yields machine precicion at Gamma_b, but requires more "radial" elements in the PML to resolve the rapid decay function
    if manuelRefinement
        pml.sigmaType = 1;  % sigmaType = 1: sigma(xi) = xi*(exp(gamma*xi)-1), sigmaType = 2: sigma(xi) = C*xi^n
        pml.dirichlet = 0;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
        pml.gamma = 2;
    else
        pml.sigmaType = 3;  % sigmaType = 1: sigma(xi) = xi*(exp(gamma*xi)-1), sigmaType = 2: sigma(xi) = C*xi^n
        pml.dirichlet = 0;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
    end
    
    switch misc.BC
        case {'SHBC','NBC'}
            k_pml = 9;
            k_P = linspace(9, 36, 3);
%             k_P = [9, 36];
%             k_P = [36, 9];
%             k_P = linspace(9, 36, 3)/5;
            k = k_P(1):0.05:k_P(end);
%             k = k_P(1):1:k_P(end);
%             k = linspace(9, 36, 5);
%             k = k_P(end);
%             k = 36;
            c_f = varCol{1}.c_f;
            f = k*c_f/(2*pi);
            f_P = k_P*c_f/(2*pi);
            omega = 2*pi*f;
            misc.omega = omega;
            omega_P = 2*pi*f_P;
            if hetmaniukCase
                misc.P_inc = -1;
                ffp.beta = pi/2;   
            end
            if msh.parm == 1
                ffp.paramPts = {[0,1,0], []};
            else
                ffp.paramPts = cell(1,12);
                ffp.paramPts{1} = [0.5,0.5,0];
            end
            for j = 1:3
                postPlot(j).plotResults = true;
                postPlot(j).printResults = true;
                postPlot(3+j).plotResults = false;
                postPlot(3+j).printResults = false;
            end

            pml.t = 0.2*varCol{1}.R;         % thickness of PML
            pml.gamma = 2.5;          % parameter for sigmaType = 1
            misc.r_a = 1.2*varCol{1}.R;
%             rom.noVecs = [8,16,24,32];
            rom.noVecs = 32;
        case {'SSBC','NNBC'}
            k_pml = 1430*2*pi/varCol{1}.c_f;
            f_P = linspace(1430, 4290, 5);
%             f_P = [1430,4290];
%             f_P = [4290,1430];
%             f_P = [4290,3146,1430];
%             f_P = [1430,3146,4290];
%             f_P = [1430,4290,3146];
%             f_P = [1430,2860,3146];
%             f_P = [1430,3146,2860];
            f_P = [1430,3146,4290];
            omega_P = 2*pi*f_P;
            f = f_P(1):12:f_P(end);
%             f = f_P(1):120:f_P(end);
%             f = f_P(1):1200:f_P(end);
%             f = linspace(1430, 4290, 9);
%             f = f_P(end);
            f = sort(unique([f_P,f]));
            omega = 2*pi*f;
            misc.omega = omega;
            k_P = omega_P/varCol{1}.c_f;
            if hetmaniukCase
                misc.P_inc = 1;
                ffp.beta = -pi/2;   
            end
            ffp.paramPts = {[0,0,0], []};
            for j = 1:3
                postPlot(j).plotResults = false;
                postPlot(j).printResults = false;
                postPlot(3+j).plotResults = true;
                postPlot(3+j).printResults = true;
            end

%             rom.noVecs = [4,8,12,20,32];
            rom.noVecs = 20;
            pml.t = 0.2*varCol{1}.R;         % thickness of PML
            pml.gamma = 2.0;          % parameter for sigmaType = 1
            misc.r_a = 1.2*varCol{1}.R;
    end
    rom.omega = omega_P;
    if ~(pml.sigmaType == 1)
        pml.gamma = 1/(k_pml*pml.t);
    end

    rom.useROM = true;

    misc.storeFullVarCol = false;
    misc.scatteringCase = 'Sweep';
    loopParameters = {'msh.M','misc.method','misc.coreMethod','misc.BC','rom.basisROM'};
    
    collectIntoTasks

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rom.useROM = false;
    for j = 1:6
        postPlot(j).noXLoopPrms   	= 1;
        switch misc.BC
            case {'SHBC','NBC'}
                postPlot(j).xname = 'varCol{1}.k';
            case {'SSBC','NNBC'}
                postPlot(j).xname = 'f';
        end
        postPlot(j).xLoopName     	= 'misc.omega';
    end
%     misc.omega = misc.omega(1);
    
    misc.scatteringCase = 'BI';
    loopParameters = {'msh.M','misc.method','misc.coreMethod','misc.BC','misc.omega'};
    collectIntoTasks
    
    %% Run paraview visualization case
    misc.omega = misc.omega(end);
    para.plotResultsInParaview = true;
    para.plotSubsets = {};
    para.plotFullDomain = 1;
    if msh.refineThetaOnly
        para.extraXiPts              = '50';  % Extra visualization points in the xi-direction per element
    else
        para.extraXiPts              = '0';  % Extra visualization points in the xi-direction per element
    end
    para.extraEtaPts             = '0';  % Extra visualization points in the eta-direction per element
    para.extraZetaPts            = '0';   % Extra visualization points in the zeta-direction per element
%     collectIntoTasks
    
    %% Run BA sweep
    misc.omega = omega;
%     misc.omega = misc.omega(1);
    para.plotResultsInParaview = 0;
    misc.method = {'BA'};
    misc.formulation = {'SL2E'};
%     misc.formulation = {'VL2E'};
    collectIntoTasks
end

function addCommands_(study)
for i_task = 1:numel(study.tasks)
    task = study.tasks(i_task).task;
    c_f = task.varCol{1}.c_f;

    if isfield(task.rom,'history')
        history = task.rom.history;
        options.xlabel = 'k';
        options.ylabel = 'residual';
        for i = 1:numel(history)
            figure(20+i)
            omega = task.rom.history(i).omega;
            residual = task.rom.history(i).residual;
            hold on
        
            omega_P = history(i).omega_P;
            J_P = history(i).J_P;
            semilogy(omega_P/c_f,1e-15*ones(size(omega_P)),'o','color','magenta','DisplayName','Interpolation points')
            options.x = omega_P.'/c_f;
            options.y = 1e-15*ones(size(omega_P)).';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_InterpolationPoints_iter' num2str(i)], options)
            for j = 1:numel(omega_P)
                text((omega_P(j)+omega_P(end)/200)/c_f,1e-15,num2str(J_P(j)),'color','magenta')
            end
            options.x = omega_P.'/c_f;
            options.y = J_P.';
            options.ylabel = 'J_P';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_noDerivatives_iter' num2str(i)], options)
    
            semilogy([omega_P(1),omega_P(end)]/c_f,task.rom.tolerance*ones(1,2),'red','DisplayName','Tolerance')
            options.x = [omega_P(1),omega_P(end)].'/c_f;
            options.y = task.rom.tolerance*ones(2,1);
            options.ylabel = 'residual';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_tolerance'], options)
    
            semilogy(omega/c_f,residual,'*','color','black','DisplayName','Residual check points')
            options.x = omega.'/c_f;
            options.y = residual.';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_residual_iter' num2str(i)], options)
    
            omega_T_new = history(i).omega_T_new;
            [~,idx] = ismember(omega,omega_T_new);
            semilogy(omega_T_new/c_f,residual(logical(idx)),'o','color','cyan','DisplayName','New interpolation point')
    
            options.x = omega_T_new.'/c_f;
            options.y = residual(logical(idx)).';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_NewInterpolationPoints_iter' num2str(i)], options)
    
            ylim([5e-16,10])
            set(gca,'yscale','log')
            legend show
            savefig([task.resultsFolder, '/', task.saveName '_adaptiveROM_iter' num2str(i)])
        end
        legend off
        if isfield(task.rom.history(end),'residualFine')
            residual = task.rom.history(end).residualFine;
            omegaFine = task.rom.history(end).omegaFine;
            semilogy(omegaFine/c_f,residual,'blue','DisplayName','Final relative residual')
        
            options.x = omegaFine.'/c_f;
            options.y = residual.';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_FinalRelativeResidual'], options)
        end
        if isfield(task.rom.history(end),'relError')
            relError = task.rom.history(end).relError;
            semilogy(omegaFine/c_f,relError,'green','DisplayName','Final relative error')
        
            options.x = omegaFine.'/c_f;
            options.y = relError.';
            options.xlabel = 'k';
            options.ylabel = 'error';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_FinalRelativeError'], options)
        end
        set(gca,'yscale','log')
        ylim([5e-16,10])
        ylabel('Relative error/residual')
        xlabel('Wavenumber')
        legend show
        savefig([task.resultsFolder, '/', task.saveName '_adaptiveROM_iter' num2str(i)])
    end
end





