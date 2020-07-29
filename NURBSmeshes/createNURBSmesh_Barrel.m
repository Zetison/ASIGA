function varCol = createNURBSmesh_Barrel(varCol, M, degree)

x_0 = [0, 0, 0];
switch varCol{1}.method
    case {'IE','IENSG','ABC'}
        varCol{1}.x_0 = x_0; % The origin of the model
        varCol{1}.A_2 = [0 1 0;
              0 0 1;
              1 0 0];
        varCol{1}.alignWithAxis = alignWithAxis;
end

if varCol{1}.boundaryMethod
    eta2 = 0.765;
    eta1 = 0.235;

    c_z = 0.99*L/2; % 30
    c_xy = 0.99*R_o; % 12
    Upsilon = sqrt(c_z^2 - c_xy^2);
    chimin = 9.8;
    chimax = 11.2;
    
    if varCol{1}.parm(1) == 1
        solid = getBarrelData(R_o,R_i,eta1,eta2,L);

        noNewXiKnots = 2^M-1;  
        noNewZetaKnots = 2^M-1;
        i_M1 = 2^M-1;
        i_M2 = 2^(M+1)-1;
        i_M3 = 2^M-1;


        solid = elevateNURBSdegree(solid,[1 1 1]*(degree-2));

        solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, noNewXiKnots) ...
                                          [insertUniform2([0 eta1], i_M1); ...
                                           insertUniform2([eta1 eta2], i_M2);
                                           insertUniform2([eta2 1], i_M3)] ...
                                          insertUniform2(solid.knots{3}, noNewZetaKnots)});



        fluid = extractOuterSurface(solid);
    else
        fluid = getBarrelData2(R_o,R_i,L);
        noNewXiKnots = 2^(M-1)-1; % 8*i_mesh
        noNewEtaKnots = noNewXiKnots;
        fluid = elevateDegreeInPatches(fluid,[1 1]*(degree-2));
        fluid = insertKnotsInPatches(fluid,noNewXiKnots,noNewEtaKnots);
    end
    varCol{1}.patchTop = getPatchTopology(fluid);
    L_gamma = L;
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGA2DMesh_new(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
else

    eta2 = 0.765;
    eta1 = 0.235;

    c_z = 15; % 30
    c_xy = c_z/2.8; % 12
    L1 = c_z/19;
    fluid = getBarrelData(R_o,R_i,eta1,eta2,L);
    intermediateLayer = getBarrelData(L1+R_o,R_i,eta1,eta2,L+2*L1);

    Upsilon = sqrt(c_z^2-c_xy^2);
    [R_a, ~, ~] = evaluateProlateCoords(c_xy,0,0,Upsilon);
    rm = zeros(N+2,1);
    for m = 1:N+2
        rm(m) = R_a + (m-1)*R_a;
    end

    fluid = embedSolidInProlateSpheroid3(fluid,intermediateLayer,c_z,c_xy,'Xaxis',x_0);

    X = evaluateNURBS(fluid, [0, eta1, 1])';
    xt = A_2*(X-x_0)';
    [~, theta_eta1, ~] = evaluateProlateCoords(xt(1),xt(2),xt(3),Upsilon);

    totArcLength1 = findArcLength(R_a,Upsilon,theta_eta1,pi);
    X = evaluateNURBS(fluid, [0, eta2, 1])';
    xt = A_2*(X-x_0)';
    [~, theta_eta3, ~] = evaluateProlateCoords(xt(1),xt(2),xt(3),Upsilon);

    X = evaluateNURBS(fluid, [0, 0.5, 1])';
    xt = A_2*(X-x_0)';
    [~, theta_eta2, ~] = evaluateProlateCoords(xt(1),xt(2),xt(3),Upsilon);

    totArcLength2 = findArcLength(R_a,Upsilon,0,theta_eta3);

    totArcLengthIntermediate1 = findArcLength(R_a,Upsilon,theta_eta2,theta_eta1);
    totArcLengthIntermediate2 = findArcLength(R_a,Upsilon,theta_eta3,theta_eta2);
    
    noNewXiKnots = ceil(0.55*c_xy*pi/2*M);  
    noNewZetaKnots = ceil(0.6364*(c_z-L/2-R_o)*M);
    i_M1 = ceil(0.55*totArcLength1*M);
    i_M2_1 = ceil(0.55*totArcLengthIntermediate1*M);
    i_M2_2 = ceil(0.55*totArcLengthIntermediate2*M);
    i_M3 = ceil(0.55*totArcLength2*M);

    findEtaKnots

    fluid = elevateNURBSdegree(fluid,[1 1 1]*(degree-2));

    fluid = insertKnotsInNURBS(fluid,{insertUniform2(fluid.knots{1}, noNewXiKnots) ...
                                      [newEtaValues1 ...
                                      insertUniform2([eta1 0.5], i_M2_1)' insertUniform2([0.5 eta2], i_M2_2)' newEtaValues2]' ...
                                      insertUniform2(fluid.knots{3}, noNewZetaKnots).^stretchingFact});


    if ~strcmp(BC, 'SHBC')
        error('This is not implemented')
    end
end
varCol{1}.nurbs = fluid;
if varCol{1}.useSolidDomain
    varCol{2}.nurbs = solid;
end
if varCol{1}.useInnerFluidDomain
    varCol{3}.nurbs = fluid_i;
end

varCol{1}.chimin = chimin;
varCol{1}.chimax = chimax;
varCol{1}.L_gamma = L_gamma;
varCol{1}.Upsilon = Upsilon;