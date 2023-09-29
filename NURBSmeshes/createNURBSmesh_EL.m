function task = createNURBSmesh_EL(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

x_0 = [0, 0, 0];
task.iem.x_0 = x_0; % The center of the model
if isfield(task.msh,'alignWithAxis')
    alignWithAxis = task.msh.alignWithAxis;
else
    alignWithAxis = 'Zaxis';
end

switch task.misc.method
    case {'IE','IENSG','ABC'}
        if ~strcmp(alignWithAxis,'Zaxis')
            warning('You may need to alter A_2 to account for this.')
        end
        task.iem.A_2 = eye(3);
        varCol{1}.alignWithAxis = alignWithAxis;
end
if isfield(task.msh, 'Xi')
    Xi = task.msh.Xi;
else
    Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
end
if isfield(varCol{1},'theta_eta')
    theta_eta = varCol{1}.theta_eta;
else
    theta_eta = pi;
end
if isfield(varCol{1},'c_z')
    c_x = varCol{1}.c_x;
    if isfield(varCol{1},'c_y')
        c_y = varCol{1}.c_y;
    else
        c_y = c_x;
    end
    c_z = varCol{1}.c_z;
    C = [c_x,c_y,c_z];
    c_min = min(C);
    [c_max,I_c_max]  = max(C);
    Upsilon = sqrt(c_max^2-c_min^2);
else
    c_x = varCol{1}.R;
    c_y = varCol{1}.R;
    c_z = varCol{1}.R;
    C = [c_x,c_y,c_z];
    c_min = min(C);
    I_c_max = 1;
    Upsilon = 0;
end
if numel(varCol) == 1
    t = eps;
else
    t = c_z - varCol{2}.R;
end
parm = task.msh.parm;
if varCol{1}.boundaryMethod
    if task.msh.refineThetaOnly
        if parm ~= 1
            error('Must have parm = 1 for pure theta refinement')
        end
        degreeVec = [2,degree,degree];
        dirs = [2,3];
    else
        degreeVec = degree*ones(1,3);
        dirs = 1:3;
    end
    varCol{1}.nurbs = getEllipsoidData('C', [c_x,c_y,c_z], 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', 0, 'Xi', Xi);
    varCol{1}.nurbs = makeUniformNURBSDegree(varCol{1}.nurbs,degreeVec(1:2));
    if numel(varCol) > 1
        varCol{2}.nurbs = getEllipsoidData('C', [c_x,c_y,c_z], 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t, 'Xi', Xi);
        varCol{2}.nurbs = makeUniformNURBSDegree(varCol{2}.nurbs,degreeVec);
    end
    c_max = max([c_x,c_y,c_z]);
    L_gamma = 2*c_max;
    if parm == 1
        refLength = c_max*pi/2;
    else
        theta_m = asec(sqrt(3));
        refLength = c_max*(pi-2*theta_m);
    end
    if isfield(varCol{1},'refinement')
        varCol{1}.nurbs = insertKnotsInNURBS(varCol{1}.nurbs,varCol{1}.refinement(M));
    end
    if numel(varCol) > 1
        if isfield(varCol{2},'refinement')
            t_fluid = c_z*2*pi/(32-pi);
            varCol{2}.nurbs = insertKnotsInNURBS(varCol{2}.nurbs,varCol{2}.refinement(M,t,t_fluid));
        end
    end
    varCol{1}.patchTop = getPatchTopology(varCol{1}.nurbs);
    if numel(varCol) > 2
        if varCol{3}.R > 0
            varCol{3}.nurbs = getEllipsoidData('C', varCol{2}.R, 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', varCol{2}.R-varCol{3}.R, 'Xi', Xi, 'theta_eta', theta_eta);
            varCol{3}.nurbs_inner = flipNURBSparametrization(subNURBS(varCol{3}.nurbs,'at',[0,0;0,0;1,0]),1);
            varCol{3}.nurbs_outer = subNURBS(varCol{2}.nurbs,'at',[0,0;0,0;1,0]);
            varCol{3}.nurbs = [varCol{3}.nurbs_inner,varCol{3}.nurbs_outer];
            
            varCol{3}.patchTop = getPatchTopology(varCol{3}.nurbs);

            varCol_varCol{3}.nurbs_inner.dimension = 1;
            varCol_varCol{3}.nurbs_inner.nurbs = varCol{3}.nurbs_inner;
            varCol_varCol{3}.nurbs_inner = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_varCol{3}.nurbs_inner)));

            varCol_varCol{3}.nurbs_outer.dimension = 1;
            varCol_varCol{3}.nurbs_outer.nurbs = varCol{3}.nurbs_outer;
            varCol_varCol{3}.nurbs_outer = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_varCol{3}.nurbs_outer)));

            varCol{3}.elemsOuter = (1:varCol_varCol{3}.nurbs_outer.noElems)+varCol_varCol{3}.nurbs_inner.noElems;
            varCol{3}.noDofsInner = varCol_varCol{3}.nurbs_inner.noDofs;
            varCol{3}.noElemsInner = varCol_varCol{3}.nurbs_inner.noElems;
        else
            varCol{3}.nurbs = subNURBS(varCol{2}.nurbs,'at',[0,0;0,0;1,0]);
            varCol{3}.nurbs = makeUniformNURBSDegree(varCol{3}.nurbs,degree);
            varCol{3}.patchTop = getPatchTopology(varCol{3}.nurbs);

            varCol_varCol{3}.nurbs.dimension = 1;
            varCol_varCol{3}.nurbs.nurbs = varCol{3}.nurbs;
            varCol_varCol{3}.nurbs = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_varCol{3}.nurbs)));
            varCol{3}.elemsOuter = 1:varCol_varCol{3}.nurbs.noElems;
            varCol{3}.noDofsInner = 0;
            varCol{3}.noElemsInner = 0;
        end
    end
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = varCol{1}.nurbs;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
else
    c_x_g = c_x; % 30
    c_y_g = c_y; % 30
    c_z_g = c_z; % 30
    C_g = [c_x_g,c_y_g,c_z_g];
    c_g_max = max(C_g);
    c_g_min = min(C_g);
    if isnan(task.misc.r_a)
        t_fluid = c_g_max*2*pi/(32-pi);
%         t_fluid = 10*R_o(1)*2*pi/(32-pi);
        c_max = c_g_max+t_fluid;
    else
        c_max = task.misc.r_a;
        if c_max < c_g_max
            error('You must set misc.r_a > c_z')
        end
        t_fluid = c_max-c_g_max;
    end
    c_min = c_g_min+t_fluid;
    Upsilon = sqrt(c_max^2-c_min^2);
    C = C+t_fluid;
    L_gamma = 2*c_z_g;

    varCol{1}.nurbs = getEllipsoidData('C', C, 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t_fluid, 'Xi', Xi, 'theta_eta', theta_eta);
    if strcmp(task.misc.method,'PML')
        C = C + task.pml.t;
    end
    c_max = max(C);
    if parm == 1
        refLength = c_max*pi/2;
    else
        theta_m = asec(sqrt(3));
        refLength = c_max*(pi-2*theta_m);
    end
    if false
        varCol{1}.nurbs = explodeNURBS(varCol{1}.nurbs,2);
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
    varCol{1}.nurbs = makeUniformNURBSDegree(varCol{1}.nurbs,degreeVec);
    
    task.misc.r_a = evaluateProlateCoords([0,0,C(end)],Upsilon);
    if numel(varCol) > 1
        varCol{2}.nurbs = getEllipsoidData('C', C_g, 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t, 'Xi', Xi, 'theta_eta', theta_eta);
        varCol{2}.nurbs = makeUniformNURBSDegree(varCol{2}.nurbs,degreeVec);
    end
    if numel(varCol) > 2
        varCol{3}.nurbs = getEllipsoidData('C', C_g - t, 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', varCol{2}.R-varCol{3}.R, 'Xi', Xi, 'theta_eta', theta_eta);
        varCol{3}.nurbs = makeUniformNURBSDegree(varCol{3}.nurbs,degreeVec);
    end
    if isfield(varCol{1},'refinement')
        varCol{1}.nurbs = insertKnotsInNURBS(varCol{1}.nurbs,varCol{1}.refinement(M));
    end
    if numel(varCol) > 1
        if isfield(varCol{2},'refinement')
            varCol{2}.nurbs = insertKnotsInNURBS(varCol{2}.nurbs,varCol{2}.refinement(M,t,t_fluid));
        end
        if isfield(task.msh,'extraSolidKnots')
            varCol{2}.nurbs = insertKnotsInNURBS(varCol{2}.nurbs,[0,0,task.msh.extraSolidKnots]);
        end
    end

    if numel(varCol) > 2
        if isfield(varCol{3},'refinement')
            varCol{3}.nurbs = insertKnotsInNURBS(varCol{3}.nurbs,varCol{3}.refinement(M));
        end
    end
end
task.msh.h_max = refLength/2^(M-1);
task.msh.dirs = dirs;
varCol{1}.refLength = refLength;
chimin = C(end)*(1-10*eps);
chimax = C(end)*(1+10*eps);
if parm == 2 && degree < 4 && (strcmp(task.misc.method,'IGA') || strcmp(task.misc.method,'C0_IGA'))
    warning(['parm=2 requires degree >= 4. Using degree=4 instead of degree=' num2str(degree)])
end
if strcmp(task.misc.method,'IE') || strcmp(task.misc.method,'IENSG')
    varCol{1}.c_z = C(end);
    varCol{1}.c_xy = C(1);
end
varCol{1}.chimax = chimax;
varCol{1}.chimin = chimin;
varCol{1}.chimax = chimax;
varCol{1}.L_gamma = L_gamma;
task.iem.Upsilon = Upsilon;
task.varCol = varCol;