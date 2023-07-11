function task = createNURBSmesh_M3(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;
if isfield(varCol{1}, 'R1')
    R1 = varCol{1}.R1;
    R2 = varCol{1}.R2;
else
    R1 = varCol{1}.R;
    R2 = varCol{1}.R;
end
if isfield(task.msh, 'Xi')
    Xi = task.msh.Xi;
else
    Xi = [0,0,0,1,1,2,2,3,3,3]/3;
end
t = varCol{1}.t;
L = varCol{1}.L;
parm = task.msh.parm;
x_0 = [0, 0, 0]; % The center of the model
alignWithAxis = 'Xaxis';
switch task.misc.method
    case {'IE','IENSG','ABC','MFS'}
        % A_2 maps the x-axis to the z-axis, the y-axis to the x-axis, and
        % the z-axis to the y-axis
        task.iem.x_0 = x_0;
        A_2 = [0 1 0;
               0 0 1;
               1 0 0];
        task.iem.A_2 = A_2;
        varCol{1}.alignWithAxis = alignWithAxis;
end

R_max = max(R1,R2);
refLength = R_max*pi/2;
if varCol{1}.boundaryMethod
    c_z = (L+R1+R2)/2;
    c_xy = sqrt(R1*R2*(L^2 - (R2 - R1)^2))/L;

    Upsilon = sqrt(c_z^2-c_xy^2);
    if isfield(varCol{1},'chimin')
        chimin = varCol{1}.chimin;
    else
        chimin = 24.4;
    end
    if isfield(varCol{1},'chimax')
        chimax = varCol{1}.chimax;
    else
        chimax = 25.7;
    end

    if task.msh.refineThetaOnly
        if parm ~= 1
            error('Must have parm = 1 for pure theta refinement')
        end
        degreeVec = [2,degree,degree];
        dirs = 2:3;
    else
        degreeVec = [degree,degree,degree];
        dirs = 1:3;
    end
    solid = getBeTSSiM3Data('R1', R1, 'R2', R2, 't', t, 'L', L, 'parm', parm, 'Xi', Xi, 'x_0', task.msh.x_0);
    solid = makeUniformNURBSDegree(solid,degreeVec);
    if numel(varCol) > 1
        varCol{2}.nurbs = solid;
    end

    varCol{1}.nurbs = subNURBS(solid, 'at', [0,0;0,0;0,1]);
    varCol{1}.patchTop = getPatchTopology(varCol{1}.nurbs);
    if numel(varCol) > 2
        varCol{3}.nurbs = subNURBS(solid, 'at', [0,0;0,0;1,0]);
        varCol{3}.patchTop = getPatchTopology(varCol{3}.nurbs);
    end
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = varCol{1}.nurbs;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
else
    t_fluid = task.misc.r_a-R2;
    if ~(strcmp(task.misc.method,'PML') || strcmp(task.misc.method,'IENSG'))
        if strcmp(task.misc.model,'MS')
            Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
            c_z = 1.2*(L+2*R1)/2; % 30
            c_xy = 1.3*R1; % 2.5, 3.75
            Upsilon = sqrt(c_z^2-c_xy^2);
        elseif strcmp(task.misc.model,'IMS')
            c_z = varCol{1}.c_z;
            c_xy = varCol{1}.c_xy;
            Upsilon = sqrt(c_z^2-c_xy^2);
        else
            s = 0.25 + 0.05*(L-R_max)/R_max;
            c_z = (L+R1+R2)/2 + s*R_max;
            c_xy = R_max + s*R_max;
            Upsilon = sqrt(c_z^2-c_xy^2);
        end
        task.misc.r_a = evaluateProlateCoords([0,0,c_z],Upsilon);
    end
    
    if task.msh.refineThetaOnly
        if parm ~= 1
            error('Must have parm = 1 for pure theta refinement')
        end
        degreeVec = [2,degree,degree];
        dirs = 2:3;
    else
        degreeVec = degree;
        dirs = 1:3;
    end
    
    if numel(varCol) > 1
        solid = getBeTSSiM3Data('R1', R1, 'R2', R2, 't', t, 'L', L, 'parm', parm, 'Xi', Xi, 'x_0', task.msh.x_0);
        solid = makeUniformNURBSDegree(solid,degreeVec);
        if parm == 1 && t(1) > t(2) % Refine around singularity
            newEtaKnots = linspace2(0,t(2)/R1(1),task.msh.extraSolidKnots);
            solid(10:15) = insertKnotsInNURBS(solid(10:15),{[],newEtaKnots,[]});
            solid(1:6) = insertKnotsInNURBS(solid(1:6),{[],1-newEtaKnots,[]});
        end

        varCol{2}.nurbs = solid;
        if isfield(task.msh,'extraSolidKnots')
            varCol{2}.nurbs = insertKnotsInNURBS(varCol{2}.nurbs,[0,0,task.msh.extraSolidKnots]);
        end
    end
    if strcmp(task.misc.method,'PML') || strcmp(task.misc.method,'IENSG')
        solid_outer = getBeTSSiM3Data('R1', R1, 'R2', R2, 't', 0, 'L', L, 'parm', parm, 'Xi', Xi, 'x_0', task.msh.x_0);
        [Gamma_a,~,c_xy,c_z] = getBeTSSiSmoothM3Data('R1', R1, 'R2', R2, 't', t_fluid, 'L', L, 'Xi', Xi, 'x_0', task.msh.x_0);
        Upsilon = sqrt(c_z^2-c_xy^2);
        chimin = 25.7;
        chimax = 27.2;
        varCol{1}.nurbs = loftNURBS({solid_outer,Gamma_a});
        varCol{1}.nurbs = makeUniformNURBSDegree(varCol{1}.nurbs,degreeVec);

        if parm == 1 && t(1) > t(2) && numel(varCol) > 1 % for compatibilty due to extra refinement around singularity in solid
            varCol{1}.nurbs(7:9) = insertKnotsInNURBS(varCol{1}.nurbs(7:9),{[],newEtaKnots,[]});
            varCol{1}.nurbs(1:3) = insertKnotsInNURBS(varCol{1}.nurbs(1:3),{[],1-newEtaKnots,[]});
        end
    end
end
task.msh.h_max = refLength/2^(M-1);
task.msh.dirs = dirs;
varCol{1}.refLength = refLength;
if parm == 2 && degree < 4
    warning(['parm=2 requires degree >= 4. Using degree=4 instead of degree=' num2str(degree)])
end
if ~strcmp(task.misc.method,'PML')
    varCol{1}.chimin = chimin;
    varCol{1}.chimax = chimax;
    varCol{1}.c_z = c_z;
    varCol{1}.c_xy = c_xy;
    task.iem.Upsilon = Upsilon;
end
varCol{1}.L_gamma = L + R1 + R2;
task.varCol = varCol;


