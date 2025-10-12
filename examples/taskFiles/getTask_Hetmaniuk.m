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
BCs = {'NNBC','SSBC','SHBC'};
if strcmp(misc.applyLoad,'pointPulsation')
    BCs = {'NBC'};
end
ffp.plotFarField = ~hetmaniukCase;
% plotFarField = true;     % If false, plots the near field instead

HetmaniukMesh = 1;
if HetmaniukMesh
    BCs = {'SHBC'};
    if 0
        misc.coreMethod = {'sub_IGA','hp_FEM'};
        msh.M = [7,6]; % 6
%         msh.M = [2,1]; % 6
    else
        misc.coreMethod = {'hp_FEM','sub_IGA'};
        msh.M = [6,7]; % 6
%         msh.M = [1,2]; % 6
        misc.coreMethod = {'sub_IGA'};
        msh.M = [7]; % 6
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

prePlot.useCamlight    = 0;        % Toggle camlight on (in plotNURBS)
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
%     collectIntoTasks
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
    connectedParameters = {{}};
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

function format_noDerivatives_results()

% Directory and file pattern
BCs = {'SHBC', 'SSBC', 'NNBC'};
iters = [2,3,5];
for j = 1:length(BCs)
    folder = 'C:\Users\jonv\results\ASIGA\Hetmaniuk2012raa\';
    filePrefix = ['Hetmaniuk2012raa_M7_PML_IGA_' BCs{j} '_basisROMDGP_noDerivatives_iter'];
    fileSuffix = '.txt';
    
    % Number of iterations/files
    numFiles = iters(j);
    
    % Initialize cell array to hold LaTeX table rows
    latexRows = cell(numFiles, 1);
    
    for i = 1:numFiles
        % Construct full file name
        filename = fullfile(folder, [filePrefix num2str(i) fileSuffix]);
    
        % Read the data (skip header row)
        data = readmatrix(filename, 'NumHeaderLines', 1);
    
        % Extract columns
        k_vals = data(:,1);
        jp_vals = data(:,2);
    
        % Convert k values to comma-separated %g format
        k_str = sprintf('%g, ', k_vals);
        k_str(end-1:end) = [];  % Remove last comma
    
        % Convert J_P values to comma-separated integer format
        jp_str = sprintf('%d, ', jp_vals);
        jp_str(end-1:end) = [];
    
        % Sum of J_P values
        total_jp = sum(jp_vals);
    
        % Format LaTeX row
        latexRows{i} = sprintf('%d\t& %s\t\t& %s\t\t& %d\t\\\\', i, k_str, jp_str, total_jp);
    end
    
    % Print LaTeX rows
    fprintf('\nGenerated LaTeX table rows for %s:\n\n', BC{1});
    for i = 1:numFiles
        fprintf('%s\n', latexRows{i});
    end
    fprintf('\n\n\n');
end

function clean_log_file()

% Input and output file paths
inputFile = 'ASIGA.log';       % <-- replace with your actual file
outputFile = 'ASIGA_clean.log';   % or use same name to overwrite

% Read the file into a cell array of lines
fid = fopen(inputFile, 'r');
lines = {};
tline = fgetl(fid);
while ischar(tline)
    lines{end+1} = tline; %#ok<AGROW>
    tline = fgetl(fid);
end
fclose(fid);

% Filter out lines that start with "Surface error = "
cleanedLines = lines(~startsWith(lines, 'Surface error = '));

% Write cleaned lines to the output file
fid = fopen(outputFile, 'w');
for i = 1:length(cleanedLines)
    fprintf(fid, '%s\n', cleanedLines{i});
end
fclose(fid);

fprintf('Cleaned file written to: %s\n', outputFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function analyze_log_file()

% File: extract_avg_timings.m


% Define cases
cases = {'SHBC', 'SSBC', 'NNBC'};

fprintf('%39s & %20s & %20s & %20s & %20s & %20s & %20s\n', 'Mesh', 'FOM avg eval time', 'FOM tot time', 'ROM avg eval time', 'ROM tot time', 'Speedup online', 'Speedup tot')

noFreqs = [541, 1201, 1201];

for j = 1:length(cases)
    BC = cases{j};
    % Read the log file
    filename = ['Hetmaniuk2012raa_M7_PML_IGA_' BC '.log'];
    fid = fopen(filename, 'r');
    logData = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    logLines = logData{1};
    
    % Initialize arrays for timing values
    totalFreqTimings = [];
    surfaceErrorTimings = [];
    totalTaskTime = NaN;
    
    % Regex patterns
    freqPattern = 'Total time spent on frequency \d+ of \d+: +([0-9.]+)';
    errorPattern = 'Calculating surface error .*?using +([0-9.]+) seconds';
    finalTaskPattern = 'Total time spent on task\s*[:=]?\s*([0-9.]+)';
    
    % Loop through lines and extract timing values
    for i = 1:length(logLines)
        line = logLines{i};
        
        % Frequency timing
        freqMatch = regexp(line, freqPattern, 'tokens');
        if ~isempty(freqMatch)
            totalFreqTimings(end + 1) = str2double(freqMatch{1}{1});
        end
        
        % Surface error timing
        errorMatch = regexp(line, errorPattern, 'tokens');
        if ~isempty(errorMatch)
            surfaceErrorTimings(end + 1) = str2double(errorMatch{1}{1});
        end
        
        % Final total task time
        finalMatch = regexp(line, finalTaskPattern, 'tokens');
        if ~isempty(finalMatch)
            totalTaskTime = str2double(finalMatch{1}{1});
        end
    end
    
    % Compute averages
    avgFreqTime = mean(totalFreqTimings);
    avgSurfaceErrorTime = mean(surfaceErrorTimings);
    FOMsweepTiming = totalTaskTime - sum(surfaceErrorTimings);
    
%     % Display results
%     fprintf('--- Timing Summary ---\n');
%     fprintf('Frequencies analyzed: %d\n', length(totalFreqTimings));
%     fprintf('Average total time per frequency: %.6f seconds\n', avgFreqTime);
%     
%     if ~isempty(surfaceErrorTimings)
%         fprintf('Surface error calculations: %d\n', length(surfaceErrorTimings));
%         fprintf('Average surface error time: %.6f seconds\n', avgSurfaceErrorTime);
%     else
%         fprintf('No "Calculating surface error" timings found.\n');
%     end
%     
%     if ~isnan(totalTaskTime)
%         fprintf('Final total time spent on task: %.6f seconds\n', totalTaskTime);
%     else
%         fprintf('Final "Total time spent on task" not found.\n');
%     end
    
    % File: extract_rom_timings.m
    
    % Read the ROM log file
    filename = ['Hetmaniuk2012raa_M7_PML_IGA_' BC '_basisROMDGP.log'];
    fid = fopen(filename, 'r');
    logData = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    logLines = logData{1};
    
    % Initialize
    surfaceErrorTimings = [];
    farfieldTimings = [];
    ROMevalTimings = [];
    adaptiveROMTime = NaN;
    ROMerrorTime = NaN;
    totalTaskTimeROM = NaN;
    
    % Regular expressions
    errorPattern = 'Calculating surface error.*?using +([0-9.]+) seconds';
    farfieldPattern = 'Computing far-field pattern.*?using +([0-9.]+) seconds';
    ROMsolPattern = 'Computing ROM solution.*?using +([0-9.]+) seconds';
    adaptiveROMPattern = 'Total time for adaptive ROM basis computations[:=]? +([0-9.]+)';
    compute_ROM_error = 'Computing final ROM residual/error .*?using +([0-9.]+) seconds';

    
    % Loop through lines and extract relevant times
    for i = 1:length(logLines)
        line = logLines{i};
    
        % Extract surface error timings
        errorMatch = regexp(line, errorPattern, 'tokens');
        if ~isempty(errorMatch)
            surfaceErrorTimings(end + 1) = str2double(errorMatch{1}{1});
        end
    
        % Extract surface error timings
        errorMatch = regexp(line, farfieldPattern, 'tokens');
        if ~isempty(errorMatch)
            farfieldTimings(end + 1) = str2double(errorMatch{1}{1});
        end
    
        % Extract surface error timings
        errorMatch = regexp(line, ROMsolPattern, 'tokens');
        if ~isempty(errorMatch)
            ROMevalTimings(end + 1) = str2double(errorMatch{1}{1});
        end
    
        % Extract total adaptive ROM basis time
        romMatch = regexp(line, adaptiveROMPattern, 'tokens');
        if ~isempty(romMatch)
            adaptiveROMTime = str2double(romMatch{1}{1});
        end
    
        % Extract total ROM error calculation time
        romErrorMatch = regexp(line, compute_ROM_error, 'tokens');
        if ~isempty(romErrorMatch)
            ROMerrorTime = str2double(romErrorMatch{1}{1});
        end
        
        % Final total task time
        finalMatch = regexp(line, finalTaskPattern, 'tokens');
        if ~isempty(finalMatch)
            totalTaskTimeROM = str2double(finalMatch{1}{1});
        end
    end
    nf = noFreqs(j);
    
    % Compute average surface error time
    avgSurfaceErrorTimeROM = mean(surfaceErrorTimings);
    avgfarfieldTimeROM = mean(farfieldTimings);
    avgROMevalTimeROM = mean(ROMevalTimings);

    ROMsweepTiming = totalTaskTimeROM - ROMerrorTime - sum(surfaceErrorTimings);
    ROMevalTiming = avgROMevalTimeROM + avgfarfieldTimeROM;
    
    % Display results
%     fprintf('--- ROM Simulation Timing Summary ---\n');
%     fprintf('Number of surface error entries: %d\n', length(surfaceErrorTimings));
%     fprintf('Average surface error time: %.6f seconds\n', avgSurfaceErrorTimeROM);
%     
%     fprintf('Total adaptive ROM basis computation time: %.6f seconds\n', adaptiveROMTime);
    speedup_online = FOMsweepTiming/nf/ROMevalTiming;
    speedup_tot = FOMsweepTiming/ROMsweepTiming;
    meshStr = sprintf('Mesh ${\\cal M}_{7,3,2}^{\\textsc{%s}}$', lower(BC));

    fprintf('%39s & %20.2f & %20.2f & %20.2f & %20.2f  & %20.2f  & %20.2f  \\\\\n', meshStr, FOMsweepTiming/nf, FOMsweepTiming, ROMevalTiming, ROMsweepTiming, speedup_online, speedup_tot)
end
% nf = 776;
% totalTaskTime_h = 79.672;
% 
% totalTaskTime = totalTaskTime_h*3600;
% ROMevalTiming = NaN;
% ROMsweepTiming_h = 42.1603;
% fprintf('%s & %20.2f & %20.2f & %20.2f & %20.2f  & %20.2f  & %20.2f  \\\\\n', 'IMS', totalTaskTime/nf, totalTaskTime_h, ROMevalTiming, ROMsweepTiming, speedup_online, speedup_tot)