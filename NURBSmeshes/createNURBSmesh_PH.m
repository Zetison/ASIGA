function task = createNURBSmesh_PH(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

L = varCol{1}.gd;
R_1 = varCol{1}.R_1;
R_2 = varCol{1}.R_2;
t = varCol{1}.t;
Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
% Xi = [0,0,0,1,1,2,2,3,3,3]/3;

x_0 = [0, 0, 0]; % The center of the model
alignWithAxis = 'Xaxis';
task.iem.x_0 = x_0;
switch task.misc.method
    case {'IE','IENSG','ABC'}
        task.iem.A_2 = [0 1 0;
                         0 0 1;
                         1 0 0];
        varCol{1}.alignWithAxis = alignWithAxis;
end
if varCol{1}.boundaryMethod
    x_1 = R_2-sqrt(R_2^2-R_1^2);
    c_z = (L+2*x_1)/2;
    c_xy = R_1; % 2.5, 3.75; 
    varCol{1}.c_z = c_z;
    varCol{1}.c_xy = c_xy;

    Upsilon = sqrt(c_z^2-c_xy^2);

    chimin = 22;
    chimax = 23;

    solid = getBeTSSiPHData('R1', R_1, 'R2', R_2, 't', t, 'gd', L,'Xi', Xi);
    solid = makeUniformNURBSDegree(solid,degree);
    solid = refineNURBSevenly(solid,(2^(M-1)-1)/(R_1*2*pi/3));

    fluid = subNURBS(solid,'at',[0 0; 0 0; 0 1]);
    varCol{1}.patchTop = getPatchTopology(fluid);
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
else
    R_max = max(R_1,R_2);
    s = 0.25 + 0.05*(L-R_max)/R_max;
    c_z = (L+R_1+R_2)/2 + s*R_max;
    c_xy = R_max + s*R_max;

    Upsilon = sqrt(c_z^2-c_xy^2);
    task.misc.r_a = evaluateProlateCoords([0,0,c_z],Upsilon);

    chimin = NaN;
    chimax = NaN;

    eta1 = 0.2;
    eta2 = 1-eta1;

    solid = getBeTSSiPHData('R1', R_1, 'R2', R_2, 't', t, 'gd', L,'Xi', Xi);
    solid = makeUniformNURBSDegree(solid,degree);

    nn = max(2^(M-1)-1,0);

    fluid = embedSolidInEllipsoid(solid,[c_z,c_xy,c_xy],alignWithAxis,x_0, 0, NaN,[0,0,0,eta1,eta1,eta2,eta2,1,1,1]);

    if M > 0
        solid(1) = insertKnotsInNURBS(solid(1),{[] [1/6 linspace2(1/3, 2/3, 8) 5/6] []});
        fluid(1) = insertKnotsInNURBS(fluid(1),{[] [eta1/2 linspace2(eta1, eta2, 8) (1+eta2)/2] []});
    end
    solid = insertKnotsInNURBS(solid,[nn+1,nn,0]);
    fluid = insertKnotsInNURBS(fluid,[nn+1,nn,nn]);
end
L_gamma = L + R_1 + R_2;

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

