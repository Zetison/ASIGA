function varCol = createNURBSmesh_M3(varCol, M, degree)

x_0 = [-L/2-(R_o1-R_o2)/2, 0, 0]; % The origin of the model
alignWithAxis = 'Xaxis';
switch varCol{1}.method
    case {'IE','IENSG','ABC','MFS'}
        % A_2 maps the x-axis to the z-axis, the y-axis to the x-axis, and
        % the z-axis to the y-axis
        A_2 = [0 1 0;
               0 0 1;
               1 0 0];
        varCol{1}.x_0 = x_0;
        varCol{1}.A_2 = A_2;
        varCol{1}.alignWithAxis = alignWithAxis;
end

if varCol{1}.boundaryMethod
    c_z = (L+R_o1+R_o2)/2;
%     c_xy = (R_o1+R_o2)/2; % 2.5, 3.75; 
    c_xy = 3*sqrt(665)*(1/20); % 2.5, 3.75; 
    
    varCol{1}.c_z = c_z;
    varCol{1}.c_xy = c_xy;

    Upsilon = sqrt(c_z^2-c_xy^2);

    chimin = 24.4;
    chimax = 25.7;

    totLength = R_o1*pi/2 + sqrt(L^2+(R_o1-R_o2)^2)+ R_o2*pi/2;
    eta1 = R_o1*pi/2/totLength;
    eta2 = (R_o1*pi/2 + sqrt(L^2+(R_o1-R_o2)^2))/totLength;
    nn = 2^(M-1)-1;

    if varCol{1}.parm(1) == 1
        solid = getModel3Data(R_o1, R_o2, t, L);
        solid = elevateNURBSdegree(solid,[1 1 1]*(degree-2));
        solid = elevateNURBSdegree(solid,[0 0 1]);
        solid = insertKnotsInNURBS(solid,{[] linspace2(eta1, eta2, 5) []});
        solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, nn) ...
                                          insertUniform2(solid.knots{2}, nn) []});
    else
        solid = getModel3Data2(R_o1, R_o2, L);
%         solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, nn) ...
%                                           insertUniform2(solid.knots{2}, nn) []});
    end



    fluid = extractSurface(solid, 'zeta', 'outer');
    varCol{1}.patchTop = getPatchTopology(fluid);
%     varCol.patchTop{1} = [ones(4,1),zeros(4,1)];
%     varCol.patchTop{1}(2,2) = NaN;
%     varCol.patchTop{1}(4,2) = NaN;
    if varCol{1}.useInnerFluidDomain
        fluid_i = extractSurface(solid, 'zeta', 'inner');
        varCol{3}.patchTop = getPatchTopology(fluid_i);
    end
else
    R_max = max(R_o1,R_o2);
    s = 0.25 + 0.05*(L-R_max)/R_max;
    c_z = (L+R_o1+R_o2)/2 + s*R_max;
    c_xy = R_max + s*R_max;

    Upsilon = sqrt(c_z^2-c_xy^2);
    varCol{1}.r_a = evaluateProlateCoords(0,0,c_z,Upsilon);

    totLength = R_o1*pi/2 + sqrt(L^2+(R_o1-R_o2)^2)+ R_o2*pi/2;
    eta1 = R_o2*pi/2/totLength;
    eta2 = (R_o2*pi/2 + sqrt(L^2+(R_o1-R_o2)^2))/totLength;

    chimin = c_z;
    chimax = 6.2;

    solid = getModel3Data(R_o1, R_o2, t, L);
    solid = elevateNURBSdegree(solid,[0 0 1]);
    solid = insertKnotsInNURBS(solid,{[] linspace2(eta1, eta2, 4) []});

    nn = 2^(M-1)-1;

    eta1m = 0.25;
    eta2m = 0.8;
%     fluid = embedSolidInEllipsoid(solid,c_z,c_xy,c_xy,alignWithAxis,x_0, 0, NaN,[eta1m,eta1m,eta2m,eta2m]);
    fluid = embedSolidInEllipsoid(solid,c_z,c_xy,c_xy,alignWithAxis,x_0, 0, NaN,[eta1m,eta1m,linspace2(eta1m,eta2m,4),eta2m,eta2m]);
%         fluid = embedSolidInEllipsoid(solid,c_z,c_xy,c_xy,alignWithAxis,x_0, 0, NaN);
    fluid = elevateNURBSdegree(fluid,[0 0 1]);
    fluid = elevateNURBSdegree(fluid,[1 1 1]*(degree-2));
%     fluid = insertKnotsInNURBS(fluid,{[] linspace2(eta1m, eta2m, 5) []});

%         fluid = insertKnotsInNURBS(fluid,{[] linspace2(eta1, eta2, 4) []});
    fluid = insertKnotsInNURBS(fluid,{insertUniform2(fluid.knots{1}, nn) ...
                                      insertUniform2(fluid.knots{2}, nn) ...
                                      insertUniform2(fluid.knots{3}, nn)});

%     fluid = insertKnotsInNURBS(fluid,{[] [eta1-eta1/(2*(nn+1)),eta2+(1-eta2)/(2*(nn+1))] []});


    solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, nn) ...
                                      insertUniform2(solid.knots{2}, nn) []});
    solid = elevateNURBSdegree(solid,[1 1 1]*(degree-2));
end
if parm == 2 && degree < 4
    warning(['parm=2 requires degree >= 4. Using degree=4 instead of degree=' num2str(degree)])
end
L_gamma = L + R_o1 + R_o2;

varCol{1}.chimin = chimin;
varCol{1}.chimax = chimax;
varCol{1}.L_gamma = L_gamma;
varCol{1}.Upsilon = Upsilon;

varCol{1}.nurbs = fluid;
if varCol{1}.useSolidDomain
    varCol{2}.nurbs = solid;
end
if varCol{1}.useInnerFluidDomain
    varCol{3}.nurbs = fluid_i;
end