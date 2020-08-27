  
switch model
    case {'SS', 'PS', 'S1', 'S1_P2', 'S3', 'S5', 'S13', 'S15', 'S35', 'IL'}
        analyticSolutionExist = true;
        isSphericalShell = true;
    otherwise
        analyticSolutionExist = false;
        isSphericalShell = false;
end
if strcmp(applyLoad,'pointPulsation') || strcmp(applyLoad,'SimpsonTorus') || strcmp(applyLoad,'radialPulsation')
    analyticSolutionExist = true;
end    
varCol{1}.analyticSolutionExist = analyticSolutionExist;
if (calculateSurfaceError || calculateVolumeError) && ~analyticSolutionExist
    error('The errors cannot be computed without an analytic solution')
end

for m = 1:numel(varCol)
    if strcmp(varCol{m}.media,'solid')
        varCol{m}.C = elasticityMatrix(varCol{m}.E,varCol{m}.nu);
    end
end

varCol{1}.isSphericalShell = isSphericalShell;

switch method
    case {'IE','ABC'}
        boundaryMethod = false;
    case {'IENSG','BEM','KDT','MFS','RT'}
        boundaryMethod = true;
    case 'BA'
        switch formulation
            case {'SL2E','SL2Etot'}
                boundaryMethod = true;
            case {'VL2E','VL2Etot'}
                boundaryMethod = false;
            otherwise
                error('Not implemented')
        end
end
varCol{1}.boundaryMethod = boundaryMethod;
for i = 1:numel(varCol)
    if mod(i,2)
        varCol{i}.boundaryMethod = boundaryMethod; % Assume only 3D elasticity is implemented
    else
        varCol{i}.boundaryMethod = false; % Assume only 3D elasticity is implemented
    end
end
    

switch scatteringCase
    case {'MS','Sweep'}
        alpha_s = alpha;
        beta_s = beta;
end

varCol{1}.coreMethod = coreMethod;
