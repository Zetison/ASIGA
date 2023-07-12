function studies = getTask_Hetmaniuk()
% This study is based on Hetmaniuk2012raa (Figure 7)
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271

counter = 1;
studies = cell(0,1);
getDefaultTaskValues
saveStudies        = false;       % save ASIGA-struct into a .mat file

%% IE simulation
hetmaniukCase = 1;

% BC = 'NNBC';
misc.model = 'Hetmaniuk2012raa';  % Spherical shell
misc.coreMethod = {'IGA','C0_IGA'};
% misc.coreMethod = {'C0_IGA'};
misc.coreMethod = {'IGA'};
% misc.coreMethod = {'hp_FEM'};
% misc.coreMethod = {'sub_IGA'};
% misc.applyLoad = 'pointPulsation';
% misc.applyLoad = 'pointCharge';
misc.applyLoad = 'planeWave';

% BCs = {'SHBC'};
% BCs = {'SSBC'};
% BCs = {'NNBC'};
BCs = {'SHBC','SSBC','NNBC'};
if strcmp(misc.applyLoad,'pointPulsation')
    BCs = {'NBC'};
end
ffp.plotFarField = ~hetmaniukCase;
% plotFarField = true;     % If false, plots the near field instead

HetmaniukMesh = 0;
if HetmaniukMesh
    BCs = {'SHBC'};
    if 0
        misc.coreMethod = {'sub_IGA','hp_FEM'};
        msh.M = [7,6]; % 6
%         msh.M = [5,4]; % 6
    else
        misc.coreMethod = {'hp_FEM','sub_IGA'};
        msh.M = [6,7]; % 6
%         msh.M = [4,5]; % 6
    end
    connectedParameters = {{'msh.M','misc.coreMethod'}};
end
if 1
    sol.solver          = 'LU';  % 'LU', 'gmres', 'cgs', 'bicgstab', 'bicgstabl', 'lsqr', 'bicg'
    sol.preconditioner  = 'diag';	% 'ilu', 'SSOR', 'CSLP', 'diag'
else
    sol.solver          = 'gmres';  % 'LU', 'gmres', 'cgs', 'bicgstab', 'bicgstabl', 'lsqr', 'bicg'
    sol.preconditioner  = 'CSLP';	% 'ilu', 'SSOR', 'CSLP', 'diag'
    sol.maxit           = 6000;     % Maximal number of iteration for iterative solver
    sol.droptol         = 1e-4;     % parameter for incomplete lu factorization with threshold and pivoting (ilutp)
    sol.beta_CSLP       = 0.1;      % parameter for the Complex Shifted Laplace Preconditioner (CSLP)
end

ffp.calculateFarFieldPattern    = true;     % Calculate far field pattern
ffp.alpha_s = 0;                            % Aspect angle of incident wave
ffp.beta_s  = -pi/2;                        % Elevation angle of incident wave
ffp.alpha   = 0;                            % Aspect angles of observation points
ffp.beta = pi/2;   
ffp.r = 1;                            % radii for near-field evaluation.
ffp.splineBasedNFPcalc = true;
warning('off','NURBS:weights')

misc.r_s = 3;

if HetmaniukMesh
    msh.parm = 2;
    msh.refineThetaOnly = 0; % only for msh.parm = 1
else
    msh.parm = 1;
    msh.refineThetaOnly = 1; % only for msh.parm = 1
    msh.autoRefine = true;
end
err.calculateSurfaceError = 1;
err.calculateVolumeError  = 0;
misc.calculateFarFieldPattern = 1;
misc.checkNURBSweightsCompatibility = false;
misc.preProcessOnly = 0;

prePlot.useCamlight    = 0;        % Toggle camlight on (in plotNURBSvec)
prePlot.plot3Dgeometry = 0;
if msh.parm == 1
    prePlot.view                = [43.314595761442661,30.547666224651493];     % Set view angle [azimuth,elevation]
else
    prePlot.view                = [77.436783356555424,9.191945753550005];     % Set view angle [azimuth,elevation]
end
prePlot.plotGeometryInfo    = false;      % Plot domain boundaries (i.e. Gamma, Gamma_a, Neumann, Dirichlet, ...)
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plotFullDomain      = 1;
% prePlot.plotControlPolygon  = 0;
% prePlot.plotSubsets         = {'xz'};
prePlot.plotSubsets         = {'Gamma'}; 
% prePlot.plotSubsets         = {}; 
% prePlot.colorFun = @(v) abs(norm2(v)-1);
if msh.parm == 1
    prePlot.resolution = [100,40,0];
    % prePlot.resolution = [400,200,0];
    % prePlot.resolution = [100,0,0];
else
    prePlot.resolution = [100,100,0];
end
% prePlot.pngResolution = '-r200';
prePlot.pngResolution = '-r800';
if msh.refineThetaOnly
    msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
else
    msh.Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
end
msh.explodeNURBS = prePlot.plot3Dgeometry && prePlot.plotFullDomain && msh.parm == 1;   % Create patches from all C^0 interfaces
if prePlot.plotFullDomain
    prePlot.format = '-png';      % Use '-png' or '-pdf' (for vector graphics)
else
    prePlot.view                = [0,0];     % Set view angle [azimuth,elevation]
    if strcmp(prePlot.plotSubsets{1},'xz')
        prePlot.format = '-pdf';      % Use '-png' or '-pdf' (for vector graphics)
    else
        prePlot.format = '-png';      % Use '-png' or '-pdf' (for vector graphics)
    end
end

misc.computeCondNumber = 0;

for i = 1:numel(BCs)
    appyCommandsAt = counter;
    misc.BC = BCs{i};
%     misc.method = {'IENSG'};
%     misc.method = {'IE'};
    misc.method = {'PML'};

    postPlot(1).xname           = 'varCol{1}.kR';
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

    postPlot(3).addCommands   	= @(study,i_study,studies) addCommands_(study);
    postPlot(6).addCommands   	= @(study,i_study,studies) addCommands_(study);
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
    msh.degree = 3;
    msh.extraSolidKnots = msh.degree;

    if strcmp(misc.method{1},'PML')
        misc.formulation = {'GSB'};
    else
        misc.formulation = {'PGC'};
    end
    if HetmaniukMesh
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
    if HetmaniukMesh
        pml.refinement = @(M) max(2^(M-4)-1,3);
    end
    if msh.refineThetaOnly
        misc.extraGP = [9-msh.degree(1),0,0];    % extra quadrature points
        ffp.extraGP = [50,0,0];    % extra quadrature points
    end
    if noDomains > 1
        if HetmaniukMesh
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
    if noDomains > 2 && HetmaniukMesh
        error('Not implemented (Case not investigated by Hetmaniuk anyways')
    end
    if ~HetmaniukMesh
        msh.M = 7; % 7
    end
    rom.basisROM = {'Pade','Taylor','DGP','Hermite','Bernstein'};  % do not put basisROMcell in loopParameters (this is done automatically)
    rom.basisROM = {'Pade','DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
    rom.basisROM = {'DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
    rom.adaptiveROM = 1;
    rom.computeROMresidualFine = ~HetmaniukMesh;
    rom.computeROMerror = ~HetmaniukMesh;
    rom.J_max = 64;
    rom.useROMconditioner = true;
    misc.symmetric = 0;

    iem.N = 50; % 9
    iem.p_ie = 5;
    iem.s_ie = 2;
    iem.IElocSup = 1;        % Toggle usage of radial shape functions in IE with local support
    
    pml.eps = 1e9*eps;      % choosing eps = eps yields machine precicion at Gamma_b, but requires more "radial" elements in the PML to resolve the rapid decay function
    if HetmaniukMesh
        pml.sigmaType = 1;  % sigmaType = 1: sigma(xi) = xi*(exp(gamma*xi)-1), sigmaType = 2: sigma(xi) = C*xi^n
        pml.dirichlet = 0;
%         pml.gamma = 2;
        pml.gamma = 2.5;
    else
        pml.sigmaType = 3;  % sigmaType = 1: sigma(xi) = xi*(exp(gamma*xi)-1), sigmaType = 2: sigma(xi) = C*xi^n
        pml.dirichlet = 1;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
    end

    c_f = varCol{1}.c_f;
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
            f = k*c_f/(2*pi);
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
            k_pml = 6;
            if false
% %             f_P = linspace(1430, 4290, 5);
%             f_P = [1430,3146,4290];
%             omega_P = 2*pi*f_P;
%             f = f_P(1):12:f_P(end);
% %             f = f_P(1):120:f_P(end);
% %             f = f_P(1):1200:f_P(end);
% %             f = linspace(1430, 4290, 9);
%             f = sort(unique([f_P,f]));
% %             f = f_P(end);
            else
                k_P = linspace(6, 18, 3);

                k = k_P(1):0.01:k_P(end);
%                 k = linspace(6, 18, 5);

                f = k*c_f/(2*pi);
            end
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
    c_f = varCol{1}.c_f;
    f_P = k_P*c_f/(2*pi);
    omega_P = 2*pi*f_P;
    rom.omega = omega_P;

%     f = [f(1),f(end)];
%     f = f(end);

    omega = 2*pi*f;
    misc.omega = omega;
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
%     for j = 1:6
%         postPlot(j).noXLoopPrms   	= 1;
%         switch misc.BC
%             case {'SHBC','NBC'}
%                 postPlot(j).xname = 'varCol{1}.k';
%             case {'SSBC','NNBC'}
%                 postPlot(j).xname = 'f';
%         end
%         postPlot(j).xLoopName     	= 'misc.omega';
%     end
%     misc.omega = misc.omega(1);
    
%     misc.scatteringCase = 'BI';
%     loopParameters = {'msh.M','misc.method','misc.coreMethod','misc.BC','misc.omega'};
    loopParameters = {'msh.M','misc.method','misc.coreMethod','misc.BC'};
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
    misc.scatteringCase = 'BI';
%     collectIntoTasks
    
    %% Run BA sweep
    misc.scatteringCase = 'Sweep';
    misc.omega = omega;
%     misc.omega = misc.omega(1);
    para.plotResultsInParaview = 0;
    misc.method = {'BA'};
    misc.formulation = {'SL2E'};
%     misc.formulation = {'VL2E'};
    collectIntoTasks
end
% Study to show the influence of refining the resolution through the
% thickness of the solid domain
for j = 1:6
    postPlot(j).plotResults = false;
    postPlot(j).printResults = false;
end
BCs = {'SSBC','NNBC'};
for i = 1:numel(BCs)
    f = 4290;
    misc.omega = 2*pi*f;
    misc.BC = BCs{i};
    switch misc.BC
        case 'SSBC'
            noDomains = 2;
        case 'NNBC'
            noDomains = 3;
    end
    varCol = setHetmaniukParameters(noDomains);
    msh.M = 7; % 7
    msh.degree = 3;
    misc.method = {'PML'};
    misc.formulation = {'GSB'};
    msh.extraSolidKnots = 0:(msh.degree+2);
    loopParameters = {'msh.M','misc.method','misc.coreMethod','misc.BC','misc.omega','msh.extraSolidKnots'};
    postPlot(7) = postPlot(1);
    postPlot(7).plotResults = true;
    postPlot(7).printResults = true;
    postPlot(7).xname           = 'msh.extraSolidKnots';
    postPlot(7).lineStyle    	= '*-';
    postPlot(7).yname        	= 'surfaceError';
    postPlot(7).noXLoopPrms   	= 1;
    postPlot(7).xLoopName       = 'msh.extraSolidKnots';
%     collectIntoTasks
    
    misc.method = {'BA'};
    misc.formulation = {'SL2E'};
    msh.extraSolidKnots = [0,(msh.degree+2)];
%     collectIntoTasks
end

function addCommands_(study)
for i_task = 1:numel(study.tasks)
    task = study.tasks(i_task).task;
    plotROMresiduals(task);
end





