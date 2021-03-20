function task = createNURBSmesh_BC(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

x_0 = [-(L+a+g2+g3)/2+a+1.5, 0, 0]; % The center of the model
alignWithAxis = 'Xaxis';
switch varCol{1}.method
    case 'IE'
        % A_2 maps the x-axis to the z-axis, the y-axis to the x-axis, and
        % the z-axis to the y-axis
        A_2 = [0 1 0;
                      0 0 1;
                      1 0 0];
        varCol{1}.x_0 = x_0;
        varCol{1}.A_2 = A_2;
    case 'IENSG'
        x_0 = [-(L+a+g2+g3)/2+a, 0, 0]; % The center of the model
        % A_2 maps the x-axis to the z-axis, the y-axis to the x-axis, and
        % the z-axis to the y-axis
        alignWithAxis = 'Xaxis';
        A_2 = [0 1 0;
                      0 0 1;
                      1 0 0];
        varCol{1}.x_0 = x_0;
        varCol{1}.A_2 = A_2;
        varCol{1}.alignWithAxis = alignWithAxis;
end

if varCol{1}.boundaryMethod
    c_z = (L+a+g2+0.5*g3)/2;
    c_xy = b;
    Upsilon = sqrt(c_z^2-c_xy^2);

    chimin = 29.3;
    chimax = 31.1;
    %%%%%%%%%%%%%

    [solid, hyp] = getBCData(parms);

    nurbsDegree = solid{1}.degree(1); % assume all degrees are equal
    solid = elevateNURBSdegree(solid,[0 0 1]);
    solid = elevateNURBSdegree(solid,[1 1 1]*(degree-nurbsDegree));

    nn = 2^(M-1)-1;

    uniqueXi = unique(solid.knots{1});
    uniqueEta = unique(solid.knots{2});
    newXiKnots = linspace2(0,uniqueXi(2), round(2*pi/3*b/hyp)-1); 
    newXiKnots = [newXiKnots,  linspace2(uniqueXi(7), 0.5, round(s/hyp)-1), 0.5];
    newXiKnots = [newXiKnots, 1-fliplr(newXiKnots(1:end-1))];

    arcLength = alpha*(g2/tan(alpha) + g2/tan((pi-alpha)/2));

    newEtaKnots = linspace2(0,uniqueEta(2), round((b - g2/tan((pi-alpha)/2) - g3*tan(alpha))/hyp)-1); 
    newEtaKnots = [newEtaKnots linspace2(uniqueEta(2),uniqueEta(3), round(g3/cos(alpha)/hyp)-1)]; 
    newEtaKnots = [newEtaKnots linspace2(uniqueEta(3),uniqueEta(4), round(arcLength/2/hyp)-1)]; 
    newEtaKnots = [newEtaKnots linspace2(uniqueEta(4),uniqueEta(5), round(arcLength/2/hyp)-1)]; 

    newEtaKnots = [newEtaKnots,  linspace2(uniqueEta(5),uniqueEta(6), round(L/hyp)-1)];
    newEtaKnots = [newEtaKnots,  linspace2(uniqueEta(6), 1, round(b*incompleteEllipticIntegral(pi/2,1-(a/b)^2)/hyp)-1)];

    insertKnots = 1;
    if insertKnots
        solid = insertKnotsInNURBS(solid,{newXiKnots newEtaKnots []});
    end
    solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, nn) ...
                                      insertUniform2(solid.knots{2}, nn) []});
    fluid = extractOuterSurface(solid);
else
    %% Set optimal parameters for radial distributions of knots (from rigid sphere analysis)

    c_z = 34; % 30
    c_xy = c_z/5; % 2.5, 3.75
    Upsilon = sqrt(c_z^2-c_xy^2);
    varCol{1}.r_a = evaluateProlateCoords(0,0,c_z,Upsilon);

    chimin = 29.3;
    chimax = 31.1;
    %%%%%%%%%%%%%

    [solid, hyp] = getBCData(parms);


    solid = elevateNURBSdegree(solid,[0 0 1]);

    Eta = solid.knots{2};
    uniqueEta = unique(Eta);

    Eta(Eta == uniqueEta(2)) = 5*uniqueEta(2);
    Eta(Eta == uniqueEta(3)) = uniqueEta(3)+0.08;
    Eta(Eta == uniqueEta(4)) = uniqueEta(4)+0.07;
    Eta(Eta == uniqueEta(5)) = uniqueEta(5)+0.06;
    Eta(Eta == uniqueEta(6)) = 0.87*uniqueEta(6);

    solid.knots{2} = Eta;
    fluid = embedSolidInEllipsoid(solid,c_z,c_xy,c_xy,alignWithAxis,x_0, -pi/2);
    fluid = elevateNURBSdegree(fluid,[0 0 1]);

    nurbsDegree = fluid{1}.degree(1); % assume all degrees are equal
    solid = elevateNURBSdegree(solid,[1 1 1]*(degree-nurbsDegree));
    fluid = elevateNURBSdegree(fluid,[1 1 1]*(degree-nurbsDegree));

    if varCol{1}.useInnerFluidDomain
        fluid_i = getBCInternalWaterData(solid,'Xaxis');
        fluid_i = elevateNURBSdegree(fluid_i,[0 0 1]);
        fluid_i = elevateNURBSdegree(fluid_i,[1 1 1]*(degree-nurbsDegree));
    end


    uniqueXi = unique(fluid.knots{1});
    uniqueEta = unique(fluid.knots{2});
    newXiKnots = linspace2(0,uniqueXi(2), round(2*pi/3*b/hyp)-1); 

    newXiKnots = [newXiKnots,  linspace2(uniqueXi(7), 0.5, round(s/hyp)-1), 0.5];
    newXiKnots = [newXiKnots, 1-fliplr(newXiKnots(1:end-1))];

    arcLength = alpha*(g2/tan(alpha) + g2/tan((pi-alpha)/2));

    newEtaKnots = [];
    newEtaKnots = [newEtaKnots uniqueEta(2)/2];  
    newEtaKnots = [newEtaKnots linspace2(uniqueEta(2),uniqueEta(3), round(1.5*g3/cos(alpha)/hyp)-1)]; 
    newEtaKnots = [newEtaKnots linspace2(uniqueEta(3),uniqueEta(4), round(1.1*arcLength/2/hyp)-1)]; 
    newEtaKnots = [newEtaKnots linspace2(uniqueEta(4),uniqueEta(5), round(1.1*arcLength/2/hyp)-1)]; 

    newEtaKnots = [newEtaKnots,  linspace2(uniqueEta(5),uniqueEta(6), round(L/hyp)-1)];
    newEtaKnots = [newEtaKnots,  linspace2(uniqueEta(6), 1, round(2*b*incompleteEllipticIntegral(pi/2,1-(a/b)^2)/hyp)-1)];

    insertKnots = 1;
    if insertKnots
        fluid = insertKnotsInNURBS(fluid,{newXiKnots newEtaKnots insertUniform2(fluid.knots{3}, 8)});
    end

    nn = 2^(M-1)-1;
    fluid = insertKnotsInNURBS(fluid,{insertUniform2(fluid.knots{1}, nn) ...
                                      insertUniform2(fluid.knots{2}, nn) ...
                                      insertUniform2(fluid.knots{3}, nn)});

    if insertKnots
        solid = insertKnotsInNURBS(solid,{newXiKnots newEtaKnots []});
    end
    solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, nn) ...
                                      insertUniform2(solid.knots{2}, nn) []});

    if varCol{1}.useInnerFluidDomain
        if insertKnots
            fluid_i = insertKnotsInNURBS(fluid_i,{newXiKnots newEtaKnots insertUniform2(fluid_i.knots{3}, 6)});
        end
        fluid_i = insertKnotsInNURBS(fluid_i,{insertUniform2(fluid_i.knots{1}, nn) ...
                                              insertUniform2(fluid_i.knots{2}, nn) []});
    end
end
L_gamma = L+a+g2+g3;


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
task.varCol = varCol;