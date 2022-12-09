function task = setAndDefineParameters(task)
switch task.misc.model
    case {'SS', 'PS', 'S1', 'S1_P2', 'S3', 'S5', 'S13', 'S15', 'S35', 'IL', 'FreeCADsphere','Shirron2006afe','Hetmaniuk2012raa','PMLstudy','Mi2021ilc'}
        analyticSolutionExist = true;
        isSphericalShell = true;
    case {'Safjan2002tdi'}
        analyticSolutionExist = true;
        isSphericalShell = (task.iem.Upsilon == 0);
    otherwise
        analyticSolutionExist = false;
        isSphericalShell = false;
end
applyLoad = task.misc.applyLoad;
if strcmp(applyLoad,'pointPulsation') || strcmp(applyLoad,'SimpsonTorus') || strcmp(applyLoad,'radialPulsation')
    analyticSolutionExist = true;
end    
if (task.err.calculateSurfaceError || task.err.calculateVolumeError) && ~analyticSolutionExist
    error('The errors cannot be computed without an analytic solution')
end

noDomains = numel(task.varCol);
for m = 1:noDomains
    if strcmp(task.varCol{m}.media,'solid')
        task.varCol{m}.C = elasticityMatrix(task.varCol{m}.E,task.varCol{m}.nu);
    end
end


boundaryMethod = isBoundaryMethod(task);

if strcmp(task.misc.scatteringCase,'MS')
    task.ffp.alpha_s = task.ffp.alpha;
    task.ffp.beta_s = task.ffp.beta;
end
task.misc.f = task.misc.omega/(2*pi);
task.omega_mean = mean([min(task.misc.omega),max(task.misc.omega)]);

for i = 1:noDomains
    switch task.varCol{i}.media
        case 'fluid'
            task.varCol{i}.k = task.misc.omega/task.varCol{i}.c_f;
            task.varCol{i}.lambda = 2*pi./task.varCol{i}.k;
            task.varCol{i}.boundaryMethod = boundaryMethod; % Assume only 3D elasticity is implemented
        case 'solid'
            task.varCol{i}.boundaryMethod = false; % Assume only 3D elasticity is implemented
    end
    if task.rom.useROM
        task.noRHSs = task.rom.noVecs;
    else
        task.noRHSs = numel(task.ffp.alpha_s);
    end
end
task.analyticSolutionExist = analyticSolutionExist;
task.isSphericalShell = isSphericalShell;
task.noDomains = noDomains;
task.boundaryMethod = boundaryMethod;
    
