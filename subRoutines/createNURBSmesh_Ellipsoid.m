function varCol = createNURBSmesh_Ellipsoid(varCol, parms, M, degree, model)

names = fieldnames(parms);
for j = 1:numel(names)
    name = names{j};
    eval([name, ' = parms.(names{j});']);
end

x_0 = [0, 0, 0];
varCol{1}.x_0 = x_0; % The origin of the model
alignWithAxis = 'Zaxis';
switch varCol{1}.method
    case {'IE','IENSG','ABC'}
        varCol{1}.A_2 = [1 0 0;
                      0 1 0;
                      0 0 1];
end

R_o = parms.R_o;
t = parms.t;
initMeshFactXi = varCol{1}.initMeshFactXi;
initMeshFactZeta = varCol{1}.initMeshFactZeta;
if varCol{1}.boundaryMethod
    switch model
        case 'EL'
            c_z = 6*R_o; % 6*R_o
            c_xy = R_o; % 2.5, 3.75
        case {'SS', 'S1', 'S1_P', 'S1_P2', 'S3', 'S5', 'SS_P','IL'}
            c_z = R_o; % 30
            c_xy = R_o; % 2.5, 3.75
    end
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
    if varCol{1}.parm(1) == 1
        solid = getEllipsoidalShellData(c_xy,c_xy,c_z,t,alignWithAxis);
        solid = elevateDegreeInPatches(solid,[0 0 1]);
    else
%         solid = getShericalShellOctPart(c_z,t);
        solid = getSphericalShellDataPatched(c_z, t); 
        solid = elevateDegreeInPatches(solid,[0 0 3]);
    end
    nurbsDegree = solid{1}.degree(1); % assume all degrees are equal
%     solid = explodeNURBS(solid,'eta');
%     solid = explodeNURBS(solid,'xi');
    degree = max(degree,nurbsDegree);
    solid = elevateDegreeInPatches(solid,[1 1 1]*(degree-nurbsDegree));
    solid = insertKnotsInPatches(solid,noNewXiKnots,noNewEtaKnots,noNewZetaKnots);
    L_gamma = 2*c_z;
    
    fluid = extractSurface(solid, 'zeta', 'outer');
    varCol{1}.patchTop = getPatchTopology(fluid);
    if varCol{1}.useInnerFluidDomain
        fluid_i = extractSurface(solid, 'zeta', 'inner');
        varCol{3}.patchTop = getPatchTopology(fluid_i);

        varCol_fluid_i.dimension = 1;
        varCol_fluid_i.nurbs = fluid_i;
        varCol_fluid_i = findDofsToRemove(generateIGA2DMesh_new(convertNURBS(varCol_fluid_i)));
        varCol{3}.elemsOuter = 1:varCol_fluid_i.noElems;
        varCol{3}.noDofsInner = 0;
        varCol{3}.noElemsInner = 0;
    end
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGA2DMesh_new(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
else
    switch model
        case 'EL'
            if 1 % if normal
                c_z_g = 6*R_o; % 30
                c_xy_g = R_o; % 2.5, 3.75
                t_fluid = c_z_g*(1.0001-1);
                c_z = c_z_g+t_fluid;
                c_xy = c_xy_g+t_fluid;
                f_arc = @(s) sqrt(c_xy_g^2*sin(s).^2+c_z_g^2*cos(s).^2);
                noNewXiKnots = 3*(2^(M-1)-1);
                noNewEtaKnots = 3*round(integral(f_arc,0,pi/2)/(c_xy_g*pi/2)*2^(M-1));
                noNewZetaKnots = 0; %2^(i_M-1);
                L_gamma = 2*c_z_g;
            else
                c_z = 0.99*(L+2*R_o)/2; % 30
                c_xy = 0.99*R_o; % 2.5, 3.75
                t_fluid = c_z/10;
            %     c_z = c_xy;
                f_arc = @(s) sqrt(c_xy^2*sin(s).^2+c_z^2*cos(s).^2);
                newXiKnots = 2^(i_M-1)-1;
                newEtaKnots = round(integral(f_arc,0,pi/2)/(c_xy*pi/2)*2^(i_M-1));
                noNewZetaKnots = round(2*t_fluid/(c_xy*pi/2)*(2^(i_M-1)+1));
                L_gamma = 2*c_z;
            end
        case {'SS', 'S1', 'S1_P', 'S1_P2', 'S3', 'S5', 'SS_P', 'IL'}
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
%                 noNewZetaKnots = max(2^(M-2),1);
            L_gamma = 2*c_z_g;
    end
    chimin = 0.99999999*c_z;
    chimax = 1.00000001*c_z;


    if varCol{1}.parm(1) == 1
        fluid = getEllipsoidalShellData(c_xy,c_xy,c_z,t_fluid,alignWithAxis);
        fluid = elevateDegreeInPatches(fluid,[0 0 1]);
    else
        fluid = getSphericalShellDataPatched(c_z, t_fluid); 
        fluid = elevateDegreeInPatches(fluid,[0 0 3]);
    end
    explodeNURBSpatches = 1;
    if explodeNURBSpatches
        fluid = explodeNURBS(fluid,'eta');
        fluid = explodeNURBS(fluid,'xi');
    end
    nurbsDegree = fluid{1}.degree(1); % assume all degrees are equal
    degree = max(degree,nurbsDegree);
    fluid = elevateDegreeInPatches(fluid,[1 1 1]*(degree-nurbsDegree));
    fluid = insertKnotsInPatches(fluid,noNewXiKnots,noNewEtaKnots,noNewZetaKnots);
    
    Upsilon = sqrt(c_z^2-c_xy^2);
    varCol{1}.r_a = evaluateProlateCoords(0,0,c_z,Upsilon);

    if strcmp(model, 'SS') || strcmp(model, 'S1') || strcmp(model, 'S3') || strcmp(model, 'S5') || strcmp(model, 'IL')
        if varCol{1}.useSolidDomain
            if varCol{1}.parm(1) == 1
                solid = getSphericalShellData(R_o, t, alignWithAxis);
                solid = elevateDegreeInPatches(solid,[0 0 1]);
            else
                solid = getSphericalShellDataPatched(R_o, t); 
                solid = elevateDegreeInPatches(solid,[0 0 3]);
            end
            if explodeNURBSpatches
                solid = explodeNURBS(solid,'eta');
                solid = explodeNURBS(solid,'xi');
            end
            nurbsDegree = solid{1}.degree(1); % assume all degrees are equal
            solid = elevateDegreeInPatches(solid,[1 1 1]*(degree-nurbsDegree));
            noNewZetaKnotsSolid = max(round(t/t_fluid*2^(M-1)),0);
%             noNewZetaKnots = 1;
            solid = insertKnotsInPatches(solid,noNewXiKnots,noNewEtaKnots,noNewZetaKnotsSolid);
        end

        if varCol{1}.useInnerFluidDomain
            if varCol{1}.parm(1) == 1
                fluid_i = getSolidSphereData(R_i, alignWithAxis);
                fluid_i = elevateDegreeInPatches(fluid_i,[0 0 1]);
            else
                fluid_i = getSphericalShellDataPatched(R_i, R_i); 
                fluid_i = elevateDegreeInPatches(fluid_i,[0 0 3]);
            end
            if explodeNURBSpatches
                fluid_i = explodeNURBS(fluid_i,'eta');
                fluid_i = explodeNURBS(fluid_i,'xi');
            end
            nurbsDegree = fluid_i{1}.degree(1); % assume all degrees are equal
            fluid_i = elevateDegreeInPatches(fluid_i,[1 1 1]*(degree-nurbsDegree));

            noNewZetaKnotsInner = max(round(R_i/t_fluid),0);
            fluid_i = insertKnotsInPatches(fluid_i,0,0,noNewZetaKnotsInner);
            fluid_i = insertKnotsInPatches(fluid_i,noNewXiKnots,noNewEtaKnots,noNewZetaKnots);
        end
    end
end