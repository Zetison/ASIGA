function varCol = createNURBSmesh_EL(varCol, M, degree)

x_0 = [0, 0, 0];
varCol{1}.x_0 = x_0; % The origin of the model
alignWithAxis = 'Zaxis';
switch varCol{1}.method
    case {'IE','IENSG','ABC'}
        varCol{1}.A_2 = [1 0 0;
                         0 1 0;
                         0 0 1];
        varCol{1}.alignWithAxis = alignWithAxis;
end

R_i = varCol{1}.R_i;
if numel(varCol) == 1
    t = R_i;
else
    t = varCol{1}.R_i - varCol{2}.R_i;
end
parm = varCol{1}.parm(1);
initMeshFactXi = varCol{1}.initMeshFactXi;
initMeshFactZeta = varCol{1}.initMeshFactZeta;
if varCol{1}.boundaryMethod
    c_z = R_i; % 30
    c_xy = R_i; % 2.5, 3.75
    % c_z = c_xy;
    if 1
        Upsilon = sqrt(c_z^2-c_xy^2);
        chimin = 0.99999999*c_z;
        chimax = 1.00000001*c_z;
    else
        Upsilon = 0;
        chimin = 0.99999999*c_xy;
        chimax = 1.00000001*c_z;
    end

    f_arc = @(s) sqrt(c_xy^2*sin(s).^2+c_z^2*cos(s).^2);
    noNewXiKnots = initMeshFactXi*2^(M-1)-1;
    noNewEtaKnots = round(integral(f_arc,0,pi/2)/(c_xy*pi/2)*(initMeshFactXi*2^(M-1)-1));
    noNewZetaKnots = max(initMeshFactZeta*2^(M-1)/8-1,0);
    solid = getEllipsoidData('C', [c_xy,c_xy,c_z], 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t);
    solid = makeUniformNURBSDegree(solid,degree);
%     solid = explodeNURBS(solid,'eta');
%     solid = explodeNURBS(solid,'xi');
    solid = insertKnotsInNURBS(solid,[noNewXiKnots,noNewEtaKnots,noNewZetaKnots]);
    L_gamma = 2*c_z;
    
    options.at = [0 0; 0 0; 0 1];
    fluid = subNURBS(solid,options);
    varCol{1}.patchTop = getPatchTopology(fluid);
    if numel(varCol) > 2
        fluid_i = extractSurface(solid, 'zeta', 'inner');
        varCol{3}.patchTop = getPatchTopology(fluid_i);

        varCol_fluid_i.dimension = 1;
        varCol_fluid_i.nurbs = fluid_i;
        varCol_fluid_i = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_fluid_i)));
        varCol{3}.elemsOuter = 1:varCol_fluid_i.noElems;
        varCol{3}.noDofsInner = 0;
        varCol{3}.noElemsInner = 0;
    end
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
else
    c_z_g = R_i; % 30
    c_xy_g = R_i; % 2.5, 3.75
    if isnan(varCol{1}.r_a)
        t_fluid = R_i(1)*2*pi/(32-pi);
%         t_fluid = 10*R_o(1)*2*pi/(32-pi);
        c_z = c_z_g+t_fluid;
        c_xy = c_xy_g+t_fluid;
    else
        c_z = varCol{1}.r_a;
        c_xy = c_z;
        t_fluid = c_z-c_z_g;
    end
    noNewXiKnots = initMeshFactXi*2^(M-1)-1;
    noNewEtaKnots = noNewXiKnots;
    noNewZetaKnots = max(initMeshFactZeta*2^(M-1)/8-1,0);
%     noNewZetaKnots = max(2^(M-2),1);
    L_gamma = 2*c_z_g;
    chimin = 0.99999999*c_z;
    chimax = 1.00000001*c_z;


    fluid = getEllipsoidData('C', [c_xy,c_xy,c_z], 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t_fluid);
    fluid = makeUniformNURBSDegree(fluid,degree);
    explodeNURBSpatches = 0;
    if explodeNURBSpatches
        fluid = explodeNURBS(fluid,'eta');
        fluid = explodeNURBS(fluid,'xi');
    end
    fluid = insertKnotsInNURBS(fluid,[noNewXiKnots,noNewEtaKnots,noNewZetaKnots]);
    
    Upsilon = sqrt(c_z^2-c_xy^2);
    varCol{1}.r_a = evaluateProlateCoords([0,0,c_z],Upsilon);

    if numel(varCol) > 1
        solid = getEllipsoidData('C', R_i, 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t);
        solid = makeUniformNURBSDegree(solid,degree);
        if explodeNURBSpatches
            solid = explodeNURBS(solid,'eta');
            solid = explodeNURBS(solid,'xi');
        end
%         noNewZetaKnotsSolid = max(2^(M-1)/2-1,0);
        noNewZetaKnotsSolid = max(round(t/t_fluid*noNewZetaKnots),0);
        solid = insertKnotsInNURBS(solid,[noNewXiKnots,noNewEtaKnots,noNewZetaKnotsSolid]);
    end

    if numel(varCol) > 2
        fluid_i = getEllipsoidData('C', varCol{2}.R_i, 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', varCol{2}.R_i-varCol{3}.R_i);
        fluid_i = makeUniformNURBSDegree(fluid_i,degree);
        if explodeNURBSpatches
            fluid_i = explodeNURBS(fluid_i,'eta');
            fluid_i = explodeNURBS(fluid_i,'xi');
        end
        noNewZetaKnotsInner = max(2^(M-1)/2-1,0);
%         noNewZetaKnotsInner = max(round(varCol{2}.R_i/t_fluid),0);
        fluid_i = insertKnotsInNURBS(fluid_i,[noNewXiKnots,noNewEtaKnots,noNewZetaKnotsInner]);
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
varCol{1}.Upsilon = Upsilon;