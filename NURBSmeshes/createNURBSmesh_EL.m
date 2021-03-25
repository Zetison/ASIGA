function task = createNURBSmesh_EL(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

x_0 = [0, 0, 0];
task.ie.x_0 = x_0; % The center of the model
alignWithAxis = 'Zaxis';
switch task.misc.method
    case {'IE','IENSG','ABC'}
        task.ie.A_2 = eye(3);
        varCol{1}.alignWithAxis = alignWithAxis;
end
if isfield(task.msh, 'Xi')
    Xi = task.msh.Xi;
else
    Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
end
if isfield(varCol{1},'c_z')
    c_x = varCol{1}.c_x;
    if ~isfield(varCol,'c_y')
        c_y = c_x;
    else
        c_y = varCol{1}.c_y;
    end
    c_z = varCol{1}.c_z;
    Upsilon = sqrt(c_z^2-c_x^2);
else
    c_x = varCol{1}.R_i;
    c_y = varCol{1}.R_i;
    c_z = varCol{1}.R_i;
    Upsilon = sqrt(c_z^2-c_x^2);
end
chimin = c_z*(1-10*eps);
chimax = c_z*(1+10*eps);
if numel(varCol) == 1
    t = eps;
else
    t = c_z - varCol{2}.R_i;
end
parm = task.msh.parm;
if varCol{1}.boundaryMethod
    solid = getEllipsoidData('C', [c_x,c_y,c_z], 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t, 'Xi', Xi);
    solid = makeUniformNURBSDegree(solid,degree);
%     solid = explodeNURBS(solid,'eta');
%     solid = explodeNURBS(solid,'xi');
    if parm == 1
        refLength = c_z*pi/2;
    else
        theta_m = asec(sqrt(3));
        refLength = c_z*(pi-2*theta_m);
    end    
    if isfield(varCol{1},'refinement')
        if numel(varCol) > 1
            t_fluid = c_z*2*pi/(32-pi);
            solid = insertKnotsInNURBS(solid,varCol{2}.refinement(M,t,t_fluid));
        else
            solid = insertKnotsInNURBS(solid,varCol{1}.refinement(M));
        end
    else
        solid = refineNURBSevenly(solid,(2^(M-1)-1)/refLength,{},0);
    end
    L_gamma = 2*c_z;
    
    options.at = [0 0; 0 0; 0 1];
    fluid = subNURBS(solid,options);
    varCol{1}.patchTop = getPatchTopology(fluid);
    if numel(varCol) > 2
        if varCol{3}.R_i > 0
            fluid_i = getEllipsoidData('C', varCol{2}.R_i, 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', varCol{2}.R_i-varCol{3}.R_i, 'Xi', Xi);
            fluid_i_inner = flipNURBSparametrization(subNURBS(fluid_i,'at',[0,0;0,0;1,0]),1);
            fluid_i_outer = subNURBS(solid,'at',[0,0;0,0;1,0]);
            fluid_i_inner = refineNURBSevenly(fluid_i_inner,(2^(M-1)-1)/(c_z*pi/2),{},0);
            fluid_i_outer = refineNURBSevenly(fluid_i_outer,(2^(M-1)-1)/(c_z*pi/2),{},0);
            fluid_i = [fluid_i_inner,fluid_i_outer];
            
            varCol{3}.patchTop = getPatchTopology(fluid_i);

            varCol_fluid_i_inner.dimension = 1;
            varCol_fluid_i_inner.nurbs = fluid_i_inner;
            varCol_fluid_i_inner = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_fluid_i_inner)));

            varCol_fluid_i_outer.dimension = 1;
            varCol_fluid_i_outer.nurbs = fluid_i_outer;
            varCol_fluid_i_outer = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_fluid_i_outer)));

            varCol{3}.elemsOuter = (1:varCol_fluid_i_outer.noElems)+varCol_fluid_i_inner.noElems;
            varCol{3}.noDofsInner = varCol_fluid_i_inner.noDofs;
            varCol{3}.noElemsInner = varCol_fluid_i_inner.noElems;
        else
            fluid_i = subNURBS(solid,'at',[0,0;0,0;1,0]);
            fluid_i = makeUniformNURBSDegree(fluid_i,degree);
            fluid_i = refineNURBSevenly(fluid_i,(2^(M-1)-1)/(c_z*pi/2),{},0);
            varCol{3}.patchTop = getPatchTopology(fluid_i);

            varCol_fluid_i.dimension = 1;
            varCol_fluid_i.nurbs = fluid_i;
            varCol_fluid_i = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_fluid_i)));
            varCol{3}.elemsOuter = 1:varCol_fluid_i.noElems;
            varCol{3}.noDofsInner = 0;
            varCol{3}.noElemsInner = 0;
        end
    end
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
else
    c_x_g = c_x; % 30
    c_y_g = c_y; % 30
    c_z_g = c_z; % 30
    Upsilon = sqrt(c_z^2-c_x^2);
    if isnan(task.misc.r_a)
        t_fluid = c_z_g*2*pi/(32-pi);
%         t_fluid = 10*R_o(1)*2*pi/(32-pi);
        c_z = c_z_g+t_fluid;
    else
        c_z = task.misc.r_a;
        t_fluid = c_z-c_z_g;
    end
    c_x = sqrt(c_z^2-Upsilon^2);
    c_y = c_x;
    L_gamma = 2*c_z_g;

    fluid = getEllipsoidData('C', [c_x,c_y,c_z], 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t_fluid, 'Xi', Xi);
    explodeNURBSpatches = 0;
    if explodeNURBSpatches
        fluid = explodeNURBS(fluid,'eta');
        fluid = explodeNURBS(fluid,'xi');
    end
    if parm == 1
        refLength = c_z*pi/2;
    else
        theta_m = asec(sqrt(3));
        refLength = c_z*(pi-2*theta_m);
    end
    fluid = makeUniformNURBSDegree(fluid,degree);
    if isfield(varCol{1},'refinement')
        fluid = insertKnotsInNURBS(fluid,varCol{1}.refinement(M));
    else
        [fluid,newKnotsIns] = refineNURBSevenly(fluid,(2^(M-1)-1)/refLength,{},0);
    end
    
    task.misc.r_a = evaluateProlateCoords([0,0,c_z],Upsilon);

    if numel(varCol) > 1
        solid = getEllipsoidData('C', [c_x_g,c_y_g,c_z_g], 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t, 'Xi', Xi);
        solid = makeUniformNURBSDegree(solid,degree);
        if explodeNURBSpatches
            solid = explodeNURBS(solid,'eta');
            solid = explodeNURBS(solid,'xi');
        end
        if isfield(varCol{2},'refinement')
            solid = insertKnotsInNURBS(solid,varCol{2}.refinement(M,t,t_fluid));
        else
            solid = refineNURBSevenly(solid,(2^(M-1)-1)/refLength,{},0,3);
            solid = insertKnotsInNURBS(solid,{newKnotsIns{1}{1},newKnotsIns{1}{2},{}});
        end
    end

    if numel(varCol) > 2
        fluid_i = getEllipsoidData('C', [c_x_g,c_y_g,c_z_g] - t, 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', varCol{2}.R_i-varCol{3}.R_i, 'Xi', Xi);
        fluid_i = makeUniformNURBSDegree(fluid_i,degree);
        if explodeNURBSpatches
            fluid_i = explodeNURBS(fluid_i,'eta');
            fluid_i = explodeNURBS(fluid_i,'xi');
        end
        if isfield(varCol{3},'refinement')
            fluid_i = insertKnotsInNURBS(fluid_i,varCol{3}.refinement(M));
        else
            fluid_i = refineNURBSevenly(fluid_i,(2^(M-1)-1)/refLength,{},0,3);
            fluid_i = insertKnotsInNURBS(fluid_i,{newKnotsIns{1}{1},newKnotsIns{1}{2},{}});
        end
    end
end
if parm == 2 && degree < 4
    warning(['parm=2 requires degree >= 4. Using degree=4 instead of degree=' num2str(degree)])
end
varCol{1}.nurbs = fluid;
if numel(varCol) > 1
    varCol{2}.nurbs = solid;
end
if numel(varCol) > 2
    varCol{3}.nurbs = fluid_i;
end

varCol{1}.chimin = chimin;
varCol{1}.chimax = chimax;
varCol{1}.L_gamma = L_gamma;
task.iem.Upsilon = Upsilon;
task.varCol = varCol;