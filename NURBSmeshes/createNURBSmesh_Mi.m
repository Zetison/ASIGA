function task = createNURBSmesh_Mi(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

x_0 = [0, 0, 0];
task.iem.x_0 = x_0; % The center of the model
alignWithAxis = 'Zaxis';
switch task.misc.method
    case {'IE','IENSG','ABC'}
        task.iem.A_2 = eye(3);
        varCol{1}.alignWithAxis = alignWithAxis;
end
if isfield(task.msh, 'Xi')
    Xi = task.msh.Xi;
else
    Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
end
if isfield(varCol{1},'c_z')
    c_x = varCol{1}.c_x;
    if ~isfield(varCol,'c_y')
        c_y = c_x;
    else
        c_y = varCol{1}.c_y;
    end
    c_z = varCol{1}.c_z;
    Upsilon = sqrt(c_z^2-c_x^2);
else
    c_x = varCol{1}.R;
    c_y = varCol{1}.R;
    c_z = varCol{1}.R;
    Upsilon = sqrt(c_z^2-c_x^2);
end
if numel(varCol) == 1
    t = eps;
else
    t = c_z - varCol{2}.R;
end
parm = task.msh.parm;

c_x_g = c_x; % 30
c_y_g = c_y; % 30
c_z_g = c_z; % 30
Upsilon = sqrt(c_z^2-c_x^2);
if isnan(task.misc.r_a)
    t_fluid = c_z_g*2*pi/(32-pi);
%         t_fluid = 10*R_o(1)*2*pi/(32-pi);
    c_z = c_z_g+t_fluid;
else
    c_z = task.misc.r_a;
    t_fluid = c_z-c_z_g;
end
c_x = sqrt(c_z^2-Upsilon^2);
c_y = c_x;
L_gamma = 2*c_z_g;

fluid = getEllipsoidData('C', [c_x,c_y,c_z], 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t_fluid, 'Xi', Xi);
if parm == 1
    refLength = c_z*pi/2;
else
    theta_m = asec(sqrt(3));
    refLength = c_z*(pi-2*theta_m);
end
if task.msh.refineThetaOnly
    if parm ~= 1
        error('Must have parm = 1 for pure theta refinement')
    end
    degreeVec = [2,degree,degree];
    dirs = 2:3;
else
    degreeVec = degree;
    dirs = 1:3;
end
fluid = makeUniformNURBSDegree(fluid,degreeVec);
fluid = insertKnotsInNURBS(fluid,{{[], [], [);

task.misc.r_a = evaluateProlateCoords([0,0,c_z],Upsilon);

if numel(varCol) > 1
    solid = getEllipsoidData('C', [c_x_g,c_y_g,c_z_g], 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', t, 'Xi', Xi);
    solid = makeUniformNURBSDegree(solid,degreeVec);
    if isfield(varCol{2},'refinement')
        solid = insertKnotsInNURBS(solid,varCol{2}.refinement(M,t,t_fluid));
    else
        solid = refineNURBSevenly(solid,(2^(M-1)-1)/refLength,{},0,3);
        solid = insertKnotsInNURBS(solid,{newKnotsIns{1}{1},newKnotsIns{1}{2},{}});
    end
end

if numel(varCol) > 2
    fluid_i = getEllipsoidData('C', [c_x_g,c_y_g,c_z_g] - t, 'alignWithAxis', alignWithAxis, 'x_0', x_0, 'parm', parm, 't', varCol{2}.R-varCol{3}.R, 'Xi', Xi);
    fluid_i = makeUniformNURBSDegree(fluid_i,degreeVec);
    if isfield(varCol{3},'refinement')
        fluid_i = insertKnotsInNURBS(fluid_i,varCol{3}.refinement(M));
    else
        fluid_i = refineNURBSevenly(fluid_i,(2^(M-1)-1)/refLength,{},0,3);
        fluid_i = insertKnotsInNURBS(fluid_i,{newKnotsIns{1}{1},newKnotsIns{1}{2},{}});
    end
end
chimin = c_z*(1-10*eps);
chimax = c_z*(1+10*eps);
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

if strcmp(task.misc.method,'IE') || strcmp(task.misc.method,'IENSG')
    varCol{1}.c_z = c_z;
    varCol{1}.c_xy = c_x;
end
varCol{1}.chimax = chimax;
varCol{1}.chimin = chimin;
varCol{1}.chimax = chimax;
varCol{1}.L_gamma = L_gamma;
task.iem.Upsilon = Upsilon;
task.varCol = varCol;