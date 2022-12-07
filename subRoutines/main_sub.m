function task = main_sub(task,loopParameters,printLog,resultsFolder)

task = controlTask(task);
task.saveName = defineFolderAndFilenames(task,loopParameters);
task.resultsFolder = resultsFolder;
task.misc.printLog = printLog;
if printLog
    fprintf('\n\nRunning the case: %s', task.saveName)
end

task = setAndDefineParameters(task);

if printLog
    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Extracting CAD data ... ')
    tic
end
task = createNURBSmesh(task);
task = collectVariables(task);
if printLog
    fprintf('using %12f seconds.', toc)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot geometry
% keyboard
if task.prePlot.plot3Dgeometry
    tic
    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Plotting geometry ... ')
    task = plotMeshAndGeometry(task);
    fprintf('using %12f seconds.', toc)
end

%% Build connectivity
tic
if printLog
    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Generating IGA mesh ... ')
end
for i = 1:task.noDomains
    task.varCol{i} = generateIGAmesh(task.varCol{i});
end

%% Find nodes to remove due to gluing and homogeneous dirichlet conditions
for i = 1:task.noDomains
    task.varCol{i} = findDofsToRemove(task.varCol{i});
end
if printLog
    fprintf('using %12f seconds.', toc)
end

%% Compute derived quantities
tic
if printLog
    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Computing derived quantities ... ')
end
task = computeDerivedQuantities(task);

%% Check NURBS compatibility
if task.misc.checkNURBSweightsCompatibility
    checkNURBSweightsCompatibility(task);
end
if ~isempty(task.prePlot.QoI)
    task.results.QoIError = integrateFunc(task, task.prePlot.QoI, task.prePlot.QoI_ref);
    fprintf('\nQoI error is %.15g\n',task.results.QoIError)
end

if printLog
    fprintf('using %12f seconds.', toc)
    fprintf('\nTotal number of elements = %d', task.totNoElems)
    fprintf('\nFinite element dofs = %d', task.FEdofs)
    fprintf('\nNumber of elements per wavelength = %.2g', min(task.varCol{1}.nepw(:)))
end

if (task.prePlot.plot3Dgeometry || task.misc.preProcessOnly) && task.prePlot.abortAfterPlotting
    return
end

%% Build stiffness matrix
t_start = tic;
if task.rom.useROM
    omega = task.misc.omega;
end

if task.rom.useROM && task.rom.adaptiveROM
    task = ASIGAassembly(task,t_start);
    task = collectMatrices(task,true,false,false);
    task = Hetmaniuk2012raa_P1(task);
else
    if task.rom.useROM
        omega_arr = task.rom.omega;
        task.V = [];
    else
        omega_arr = task.misc.omega;
    end
    if task.rom.useROM && strcmp(task.misc.scatteringCase,'Sweep')
        task.U_sweep = cell(1,numel(omega_arr));
    end
    if ~(strcmp(task.misc.method,'RT') || strcmp(task.misc.method,'KDT'))
        for i_o = 1:numel(omega_arr)
            omega_i = omega_arr(i_o);
            t_freq = tic;
            task.misc.omega = omega_i;
            task = ASIGAassembly(task,t_start);
            task = collectMatrices(task);
            task = ASIGAsolve(task,omega_arr,i_o);
        
            if task.misc.printLog
                fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Total time building system ... ')
                fprintf(' was  %12f seconds.', task.timeBuildSystem+task.timeCollectingMatrices)
            end
            if printLog && strcmp(task.misc.scatteringCase,'Sweep') && numel(omega_arr) > 1
                fprintf('\nTotal time spent on frequency %d of %d: %12f\n', i_o, numel(omega_arr), toc(t_freq))  
            end
        end
        
        if ~task.rom.useROM
            task = postProcessSolution(task);
            if task.err.calculateSurfaceError || task.err.calculateVolumeError
                task = calculateErrors(task, printLog, task.misc.stringShift);
            end
        end
    else
        task.varCol{1}.U = [];
    end
end
if task.rom.useROM
    task.misc.omega = omega;
    task = computeROMsolution(task,printLog);
end
%% Compute scattered pressure   
if task.ffp.calculateFarFieldPattern && ~task.rom.useROM
    task = getAnalyticSolutions(task);
    task = calculateTS(task,printLog,task.misc.stringShift);
    task = computeDerivedFFPquantities(task,task.ffp.p_h);
end
%% Calculate errors (if analyticSolutionExist) and plot result in Paraview
task = ASIGAparaview(task);

task.varCol{1}.tot_time = toc(t_start);

if printLog
    fprintf('\nTotal time spent on task: %12f', task.varCol{1}.tot_time)  
end

if ~task.misc.storeFullVarCol
    task.varCol = rmfields(task.varCol,getAddedFields());
end






