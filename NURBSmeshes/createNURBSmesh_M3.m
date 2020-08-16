function varCol = createNURBSmesh_M3(varCol, M, degree)
R1 = varCol{1}.R1;
R2 = varCol{1}.R2;
t = varCol{1}.t;
L = varCol{1}.L;
parm = varCol{1}.parm;
x_0 = [-L/2-(R2-R1)/2, 0, 0]; % The origin of the model
alignWithAxis = 'Xaxis';
switch varCol{1}.method
    case {'IE','IENSG','ABC','MFS'}
        % A_2 maps the x-axis to the z-axis, the y-axis to the x-axis, and
        % the z-axis to the y-axis
        varCol{1}.x_0 = x_0;
        varCol{1}.A_2 = [0 1 0;
                         0 0 1;
                         1 0 0];
        varCol{1}.alignWithAxis = alignWithAxis;
end

if varCol{1}.boundaryMethod
    c_z = (L+R1+R2)/2;
    c_xy = sqrt(R1*R2*(L^2 - (R2 - R1)^2))/L;
    varCol{1}.c_z = c_z;
    varCol{1}.c_xy = c_xy;

    Upsilon = sqrt(c_z^2-c_xy^2);

    chimin = 24.4;
    chimax = 25.7;

    solid = getBeTSSiM3Data('R1', R1, 'R2', R2, 't', t, 'L', L, 'parm', parm);
    solid = makeUniformNURBSDegree(solid,degree);
    Imap{1} = [R2*pi/2,R1*pi/2,(R2-t)*pi/2,(R1-t)*pi/2];
    solid = refineNURBSevenly(solid,(2^(M-1)-1)/(R2*pi/2),Imap);

    fluid = subNURBS(solid, 'at', [0,0;0,0;0,1]);
    varCol{1}.patchTop = getPatchTopology(fluid);
    if numel(varCol) > 2
        fluid_i = subNURBS(solid, 'at', [0,0;0,0;1,0]);
        varCol{3}.patchTop = getPatchTopology(fluid_i);
    end
else
    error('The mesh is not fully implemented')
    R_max = max(R1,R2);
    s = 0.25 + 0.05*(L-R_max)/R_max;
    c_z = (L+R1+R2)/2 + s*R_max;
    c_xy = R_max + s*R_max;

    Upsilon = sqrt(c_z^2-c_xy^2);
    varCol{1}.r_a = evaluateProlateCoords([0,0,c_z],Upsilon);

    totLength = R1*pi/2 + sqrt(L^2+(R1-R2)^2)+ R2*pi/2;
    eta1 = R2*pi/2/totLength;
    eta2 = (R2*pi/2 + sqrt(L^2+(R1-R2)^2))/totLength;

    chimin = c_z;
    chimax = 6.2;
    
    solid = getBeTSSiM3Data('R1', R1, 'R2', R2, 't', t, 'L', L, 'parm', parm);
    solid = makeUniformNURBSDegree(solid,degree);
    Imap{1} = [R2*pi/2,R1*pi/2,(R2-t)*pi/2,(R1-t)*pi/2];
    solid = refineNURBSevenly(solid,(2^(M-1)-1)/(R2*pi/2),Imap);

    nn = 2^(M-1)-1;

    eta1m = 0.25;
    eta2m = 0.8;
%     fluid = embedSolidInEllipsoid(solid,[c_z,c_xy,c_xy],alignWithAxis,x_0, 0, NaN,[0,0,0,eta1m,eta1m,linspace2(eta1m,eta2m,4),eta2m,eta2m,1,1,1]);
%     fluid = embedSolidInEllipsoid(solid,[c_z,c_xy,c_xy],alignWithAxis,x_0, 0, NaN);
    temp = getFreeNURBS(solid);
%     fluid = refineNURBSevenly(fluid,2.2);
%     fluid = insertKnotsInNURBS(fluid,(2^(M-1)-1)*[1,1]);
    nurbsCol = collectConnectedNURBS(temp);
    solid_o = nurbsCol{1};
    
    ellipsoid = getEllipsoidData('C',[c_z,c_xy,c_xy],'parm',parm,'alignWithAxis',alignWithAxis,'x_0',x_0);
    
    fluid = loftNURBS(solid_o,ellipsoid);
    fluid = makeUniformNURBSDegree(fluid,degree);
end
if parm == 2 && degree < 4
    warning(['parm=2 requires degree >= 4. Using degree=4 instead of degree=' num2str(degree)])
end
varCol{1}.chimin = chimin;
varCol{1}.chimax = chimax;
varCol{1}.L_gamma = L + R1 + R2;
varCol{1}.Upsilon = Upsilon;

varCol{1}.nurbs = fluid;
if numel(varCol) > 1
    varCol{2}.nurbs = solid;
end
if numel(varCol) > 2
    varCol{3}.nurbs = fluid_i;
end