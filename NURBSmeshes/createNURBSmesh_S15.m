function task = createNURBSmesh_S15(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

names = fieldnames(parms);
for j = 1:numel(names)
    name = names{j};
    eval([name, ' = parms.(names{j});']);
end

x_0 = [0, 0, 0];
varCol{1}.x_0 = x_0; % The center of the model
alignWithAxis = 'Zaxis';
switch varCol{1}.method
    case {'IE','IENSG','ABC'}
        varCol{1}.A_2 = [1 0 0;
                      0 1 0;
                      0 0 1];
        varCol{1}.alignWithAxis = alignWithAxis;
end

R_o = parms.R_o;
t = parms.t;
parm = varCol{1}.parm(1);
initMeshFactXi = varCol{1}.initMeshFactXi;
initMeshFactZeta = varCol{1}.initMeshFactZeta;
if varCol{1}.boundaryMethod
    noNewXiKnots = initMeshFactXi*2^(M-1)-1;
    noNewEtaKnots = noNewXiKnots;
    noNewZetaKnots = max(initMeshFactZeta*2^(M-1)/8-1,0);
    solid = getEllipsoidalData({'R', R_o(1), 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t(1)});
    solid = makeUniformNURBSDegree(solid);
    solid2 = getEllipsoidalData({'R', R_o(2), 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t(2)});
    solid2 = makeUniformNURBSDegree(solid2);
    
    nurbsDegree = solid{1}.degree(1); % assume all degrees are equal
%     solid = explodeNURBS(solid,'eta');
%     solid = explodeNURBS(solid,'xi');
    degree = max(degree,nurbsDegree);
    solid = elevateDegreeInPatches(solid,[1 1 1]*(degree-nurbsDegree));
    solid2 = elevateDegreeInPatches(solid2,[1 1 1]*(degree-nurbsDegree));
    solid = insertKnotsInPatches(solid,noNewXiKnots,noNewEtaKnots,noNewZetaKnots);
    solid2 = insertKnotsInPatches(solid2,noNewXiKnots,noNewEtaKnots,noNewZetaKnots);
    L_gamma = 2*R_o(1);

    fluid = extractSurface(solid, 'zeta', 'outer');
    varCol{1}.patchTop = getPatchTopology(fluid);
    if varCol{1}.useInnerFluidDomain
        fluid_i_inner = flipNURBSparametrization(extractSurface(solid2, 'zeta', 'outer'),'xi');
        fluid_i_outer = extractSurface(solid, 'zeta', 'inner'); 
        fluid_i = uniteNURBS({fluid_i_inner, fluid_i_outer});
        varCol{3}.patchTop = getPatchTopology(fluid_i);

        varCol_fluid_i_inner.dimension = 1;
        varCol_fluid_i_inner.nurbs = fluid_i_inner;
        varCol_fluid_i_inner = findDofsToRemove(generateIGA2DMesh_new(convertNURBS(varCol_fluid_i_inner)));
        
        varCol_fluid_i_outer.dimension = 1;
        varCol_fluid_i_outer.nurbs = fluid_i_outer;
        varCol_fluid_i_outer = findDofsToRemove(generateIGA2DMesh_new(convertNURBS(varCol_fluid_i_outer)));
        
        varCol{3}.elemsOuter = (1:varCol_fluid_i_outer.noElems)+varCol_fluid_i_inner.noElems;
        varCol{3}.noDofsInner = varCol_fluid_i_inner.noDofs;
        varCol{3}.noElemsInner = varCol_fluid_i_inner.noElems;
    end
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGA2DMesh_new(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
else
    c_z_g = R_o(1); % 30
    c_xy_g = R_o(1); % 2.5, 3.75
    if isfield(varCol{1},'r_a')
        c_z = varCol{1}.r_a;
        c_xy = c_z;
        t_fluid = c_z-c_z_g;
    else
        t_fluid = R_o(1)*2*pi/(32-pi);
%                 t_fluid = 10*R_o(1)*2*pi/(32-pi);
        c_z = c_z_g+t_fluid;
        c_xy = c_xy_g+t_fluid;
    end
    noNewXiKnots = initMeshFactXi*2^(M-1)-1;
    noNewEtaKnots = noNewXiKnots;
    noNewZetaKnots = max(initMeshFactZeta*2^(M-1)/8-1,0);
    L_gamma = 2*c_z_g;


    fluid = getEllipsoidalData({'R', [c_xy,c_xy,c_z], 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t_fluid});
    fluid = makeUniformNURBSDegree(fluid);
    
    nurbsDegree = fluid{1}.degree(1); % assume all degrees are equal
    degree = max(degree,nurbsDegree);
    fluid = elevateDegreeInPatches(fluid,[1 1 1]*(degree-nurbsDegree));
    fluid = insertKnotsInPatches(fluid,noNewXiKnots,noNewEtaKnots,noNewZetaKnots);
    
    Upsilon = sqrt(c_z^2-c_xy^2);
    varCol{1}.r_a = evaluateProlateCoords(0,0,c_z,Upsilon);
    varCol{1}.Upsilon = Upsilon;
    
    if numel(varCol) > 1
        solid = getEllipsoidalData({'R', R_o(1), 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t(1)});
        solid = makeUniformNURBSDegree(solid);
        
        nurbsDegree = solid{1}.degree(1); % assume all degrees are equal
        solid = elevateDegreeInPatches(solid,[1 1 1]*(degree-nurbsDegree));
        noNewZetaKnotsSolid = max(round(t(1)/t_fluid*2^(M-1)),0);
        
        solid = insertKnotsInPatches(solid,noNewXiKnots,noNewEtaKnots,noNewZetaKnotsSolid);
    end

    if numel(varCol) > 2
        t_ifluid = R_i(1)-R_o(2);
        solid = getEllipsoidalData({'R', R_i(1), 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t_ifluid});
        solid = makeUniformNURBSDegree(solid);
        nurbsDegree = fluid_i{1}.degree(1); % assume all degrees are equal
        fluid_i = elevateDegreeInPatches(fluid_i,[1 1 1]*(degree-nurbsDegree));

        noNewZetaKnotsInner = max(round(t_ifluid/t_fluid),0);
        fluid_i = insertKnotsInPatches(fluid_i,0,0,noNewZetaKnotsInner);
        fluid_i = insertKnotsInPatches(fluid_i,noNewXiKnots,noNewEtaKnots,noNewZetaKnots);
    end
end
varCol{1}.nurbs = fluid;
if numel(varCol) > 1
    varCol{2}.nurbs = solid;
end
if numel(varCol) > 2
    varCol{3}.nurbs = fluid_i;
end

varCol{1}.L_gamma = L_gamma;
task.varCol = varCol;