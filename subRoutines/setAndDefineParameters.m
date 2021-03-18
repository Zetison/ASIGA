function task = setAndDefineParameters(task)
switch task.misc.model
    case {'SS', 'PS', 'S1', 'S1_P2', 'S3', 'S5', 'S13', 'S15', 'S35', 'IL', 'FreeCADsphere'}
        analyticSolutionExist = true;
        isSphericalShell = true;
    case {'Safjan2002tdi'}
        analyticSolutionExist = true;
        isSphericalShell = false;
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

for m = 1:numel(task.varCol)
    if strcmp(task.varCol{m}.media,'solid')
        task.varCol{m}.C = elasticityMatrix(task.varCol{m}.E,task.varCol{m}.nu);
    end
end


switch task.misc.method
    case {'IE','ABC','PML'}
        boundaryMethod = false;
    case {'IENSG','BEM','KDT','MFS','RT'}
        boundaryMethod = true;
    case 'BA'
        switch task.misc.formulation
            case {'SL2E','SL2Etot'}
                boundaryMethod = true;
            case {'VL2E','VL2Etot'}
                boundaryMethod = false;
            otherwise
                error('Not implemented')
        end
end

task.analyticSolutionExist = analyticSolutionExist;
task.isSphericalShell = isSphericalShell;

for i = 1:numel(task.varCol)
    if mod(i,2)
        task.varCol{i}.boundaryMethod = boundaryMethod; % Assume only 3D elasticity is implemented
    else
        task.varCol{i}.boundaryMethod = false; % Assume only 3D elasticity is implemented
    end
    if task.rom.useROM
        task.noRHSs = max(task.rom.noVecsArr);
    else
        task.noRHSs = numel(task.ffp.alpha_s);
    end
end
    
