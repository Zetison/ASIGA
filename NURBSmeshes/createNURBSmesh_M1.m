function task = createNURBSmesh_M1(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;


R = varCol{1}.R;
t = varCol{1}.t;
L = varCol{1}.L;

x_0 = [0, 0, 0]; % The center of the model
alignWithAxis = 'Xaxis';
switch varCol{1}.method
    case {'IE','IENSG','MFS'}
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
    c_z = 0.98*(L+R)/2;
    c_xy = 0.99*R/2; % 2.5, 3.75; 
    Upsilon = sqrt(c_z^2-c_xy^2);
    varCol{1}.c_z = c_z;
    varCol{1}.c_xy = c_xy;

    chimin = 21.07;
    chimax = 23.1;
    
    fluid = getBeTSSiM1Data('R', R, 'L', L, 't', t);
    fluid = makeUniformNURBSDegree(fluid,degree);
    Imap{1} = [R*pi/2 4.242171326235288];
    Imap{2} = [R*pi/4 2.121085663117643];
    fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/(R*pi/2),Imap);
    varCol{1}.patchTop = getPatchTopology(fluid);
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
end
varCol{1}.nurbs = fluid;
if numel(varCol) > 1
    varCol{2}.nurbs = solid;
end
if numel(varCol) > 2
    varCol{3}.nurbs = fluid_i;
end

varCol{1}.L_gamma = L + R;
varCol{1}.Upsilon = Upsilon;
varCol{1}.chimin = chimin;
varCol{1}.chimax = chimax;
task.varCol = varCol;

