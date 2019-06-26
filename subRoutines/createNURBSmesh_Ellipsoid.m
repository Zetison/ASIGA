function [varCol, fluid, solid, fluid_i] = createNURBSmesh_Ellipsoid(varCol, parms, M, degree, model)

solid = NaN;
fluid_i = NaN;

names = fieldnames(parms);
for j = 1:numel(names)
    name = names{j};
    eval([name, ' = parms.(names{j});']);
end

x_0 = [0, 0, 0];
varCol.x_0 = x_0; % The origin of the model
alignWithAxis = 'Zaxis';
switch varCol.method
    case {'IE','IENSG','ABC'}
        varCol.A_2 = [1 0 0;
                      0 1 0;
                      0 0 1];
end

R_o = parms.R_o;
t = parms.t;
initMeshFactXi = varCol.initMeshFactXi;
initMeshFactZeta = varCol.initMeshFactZeta;
if varCol.boundaryMethod
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

    t = c_z/10;

    f_arc = @(s) sqrt(c_xy^2*sin(s).^2+c_z^2*cos(s).^2);
    noNewXiKnots = initMeshFactXi*2^(M-1)-1;
    noNewEtaKnots = round(integral(f_arc,0,pi/2)/(c_xy*pi/2)*(initMeshFactXi*2^(M-1)-1));
    noNewZetaKnots = initMeshFactZeta*2^(M-1);
    if varCol.parm(1) == 1
        solid = getEllipsoidalShellData(c_xy,c_xy,c_z,t,alignWithAxis);
%         varCol.patchTop{1} = [1, 0;
%                               NaN, NaN;
%                               1, 0;
%                               NaN, NaN];
        nurbsDegree = solid.degree(1); % assume all degrees are equal
    else
        solid = getSphericalShellDataPatched(c_z, t); 
        solid = elevateDegreeInPatches(solid,[0 0 3]);
%         varCol.patchTop = cell(6,1);
%         varCol.patchTop{1} = [2 0;
%                               6 0;
%                               4 0;
%                               5 0];
%         varCol.patchTop{2} = [3 0;
%                               6 3;
%                               1 0;
%                               5 1];
%         varCol.patchTop{3} = [4 0;
%                               6 2;
%                               2 0;
%                               5 2];
%         varCol.patchTop{4} = [1 0;
%                               6 1;
%                               3 0;
%                               5 3];
%         varCol.patchTop{5} = [2 3;
%                               1 0;
%                               4 1;
%                               3 2];
%         varCol.patchTop{6} = [2 1;
%                               3 2;
%                               4 3;
%                               1 0];
        nurbsDegree = solid{1}.degree(1); % assume all degrees are equal
    end
%     solid = explodeNURBS(solid,'eta');
%     solid = explodeNURBS(solid,'xi');
    degree = max(degree,nurbsDegree);
    solid = elevateDegreeInPatches(solid,[1 1 1]*(degree-nurbsDegree));
    solid = insertKnotsInPatches(solid,noNewXiKnots,noNewEtaKnots,noNewZetaKnots);
    L_gamma = 2*c_z;

    fluid = extractOuterSurface(solid);
    varCol.patchTop = getPatchTopology(fluid);
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
        case {'SS', 'S1', 'S1_P', 'S1_P2', 'S2', 'S3', 'SS_P', 'IL'}
            c_z_g = R_o(1); % 30
            c_xy_g = R_o(1); % 2.5, 3.75
            if isfield(varCol,'r_a')
                c_z = varCol.r_a;
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


    if varCol.parm(1) == 1
        fluid = getEllipsoidalShellData(c_xy,c_xy,c_z,t_fluid,alignWithAxis);
        fluid = elevateDegreeInPatches(fluid,[0 0 1]);
    else
        fluid = getSphericalShellDataPatched(c_z, t_fluid); 
        fluid = elevateDegreeInPatches(fluid,[0 0 3]);
    end
%     fluid = explodeNURBS(fluid,'eta');
%     fluid = explodeNURBS(fluid,'xi');
    nurbsDegree = fluid{1}.degree(1); % assume all degrees are equal
    degree = max(degree,nurbsDegree);
    fluid = elevateDegreeInPatches(fluid,[1 1 1]*(degree-nurbsDegree));
    fluid = insertKnotsInPatches(fluid,noNewXiKnots,noNewEtaKnots,noNewZetaKnots);
    
    Upsilon = sqrt(c_z^2-c_xy^2);
    varCol.r_a = evaluateProlateCoords(0,0,c_z,Upsilon);

    if strcmp(model, 'SS') || strcmp(model, 'S1') || strcmp(model, 'S2') || strcmp(model, 'S3') || strcmp(model, 'IL')
        if varCol.useSolidDomain
            solid = getSphericalShellData(R_o, t, alignWithAxis);
            solid = elevateNURBSdegree(solid,[0 0 1]);
            solid = elevateNURBSdegree(solid,[1 1 1]*(degree-nurbsDegree));
            noNewZetaKnots = max(round(t/t_fluid)*2^(M-1),0);
            solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, noNewXiKnots) ...
                                              insertUniform2(solid.knots{2}, noNewEtaKnots) ...
                                              insertUniform2(solid.knots{3}, noNewZetaKnots)});
        end

        if varCol.useInnerFluidDomain
            fluid_i = getSolidSphereData(R_i, alignWithAxis);
            fluid_i = elevateNURBSdegree(fluid_i,[0 0 1]);
            fluid_i = elevateNURBSdegree(fluid_i,[1 1 1]*(degree-nurbsDegree));

            noNewZetaKnots = max(2^(M-1)/2-1,0);
            fluid_i = insertKnotsInNURBS(fluid_i,{insertUniform2(fluid_i.knots{1}, noNewXiKnots) ...
                                              insertUniform2(fluid_i.knots{2}, noNewEtaKnots) ...
                                              insertUniform2(fluid_i.knots{3}, noNewZetaKnots)});
        end
    end
end

varCol.chimin = chimin;
varCol.chimax = chimax;
varCol.L_gamma = L_gamma;
varCol.Upsilon = Upsilon;