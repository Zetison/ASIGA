function task = createNURBSmesh_Barrel(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

x_0 = [0, 0, 0];
switch varCol{1}.method
    case {'IE','IENSG','ABC'}
        varCol{1}.x_0 = x_0; % The origin of the model
        varCol{1}.A_2 = [0 1 0;
              0 0 1;
              1 0 0];
        varCol{1}.alignWithAxis = alignWithAxis;
end
Xi = [0,0,0,1,1,2,2,3,3,3]/3;
R = varCol{1}.R;
t = varCol{1}.t;
L = varCol{1}.L;
parm = varCol{1}.parm;
if varCol{1}.boundaryMethod
    c_z = L/2; % 30
    c_xy = R; % 12
    Upsilon = sqrt(c_z^2 - c_xy^2);
    chimin = 9.8;
    chimax = 11.2;

    if numel(varCol) == 1
        fluid = getBarrelData('R', R, 't', t, 'parm', parm, 'L', L, 'd_p', 2, 'Xi', Xi);
        fluid = makeUniformNURBSDegree(fluid,degree);
        
        if parm == 1
            Imap{1} = [];
            fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/(R*2*pi/3),Imap);
        else
            Imap{1} = [R*pi/2, R*1.414057108745095];
            fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/(R*pi/2),Imap);
        end
    else
        solid = getBarrelData('R', R, 't', t, 'parm', parm, 'L', L, 'Xi', Xi);
        solid = makeUniformNURBSDegree(solid,degree);
        fluid = getFreeNURBS(solid);
        if parm == 1
            Imap{1} = [R*2*pi/3, (R-t)*2*pi/3, R*1.272651397870586, R*1.155553578414605];
            fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/(R*2*pi/3),Imap);
        else
            Imap{1} = [R*pi/2, (R-t)*pi/2, R*1.272651397870586, R*1.155553578414605];
            fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/(R*pi/2),Imap);
        end
        nurbsCol = collectConnectedNURBS(fluid);
        fluid = nurbsCol{1};
        fluid_i = nurbsCol{2};
    end
    
    varCol{1}.patchTop = getPatchTopology(fluid);
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
else
    error('This is not implemented')

    eta2 = 0.765;
    eta1 = 0.235;

    c_z = 15; % 30
    c_xy = c_z/2.8; % 12
    L1 = c_z/19;
    fluid = getBarrelData(R,R_i,eta1,eta2,L);
    intermediateLayer = getBarrelData(L1+R,R_i,eta1,eta2,L+2*L1);

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
    noNewZetaKnots = ceil(0.6364*(c_z-L/2-R)*M);
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
varCol{1}.L_gamma = L;
varCol{1}.Upsilon = Upsilon;
task.varCol = varCol;