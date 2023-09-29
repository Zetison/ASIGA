function task = createNURBSmesh_Barrel(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

x_0 = [0, 0, 0];
switch task.misc.method
    case {'IE','IENSG','ABC'}
        task.iem.x_0 = x_0; % The center of the model
        task.iem.A_2 = [0 1 0;
                      0 0 1;
                      1 0 0];
        varCol{1}.alignWithAxis = 'Xaxis';
end

if isfield(varCol{1},'Xi')
    Xi = varCol{1}.Xi;
else
    Xi = [0,0,0,1,1,2,2,3,3,3]/3;
end
R = varCol{1}.R;
t = varCol{1}.t;
L = varCol{1}.L;
parm = task.msh.parm;
if varCol{1}.boundaryMethod
    c_z = L/2; % 30
    c_xy = R; % 12
    Upsilon = sqrt(c_z^2 - c_xy^2);
    chimin = 1.5;
    chimax = 2.1;
    d_p = 2;
    if task.msh.refineThetaOnly
        dirs = 2:d_p;
    else
        dirs = 1:d_p;
    end

    if numel(varCol) == 1
        fluid = getBarrelData('R', R, 't', t, 'parm', parm, 'L', L, 'd_p', 2, 'Xi', Xi);
        fluid = makeUniformNURBSDegree(fluid,degree);
        
        if parm == 1
            Imap{1} = [];
            fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/(R*2*pi/3),Imap,0,dirs);
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
            fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/(R*2*pi/3),Imap,0,dirs);
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
    d_p = 3;
    if task.msh.refineThetaOnly
        dirs = 2:d_p;
    else
        dirs = 1:d_p;
    end
    if numel(varCol) > 1
        error('not implemented')
    end
    if numel(Xi) == 12
        xiAngle = pi/2;
    else
        xiAngle = 2*pi/3;
    end
    t_fluid = task.misc.r_a - R;
    t_pml = task.pml.t;
    c_z = L/2+t_fluid;
    c_xy = R+t_fluid; % 12
    Upsilon = sqrt(c_z^2 - c_xy^2);
    chimin = 1.8;
    chimax = 2.7;
    fluid = getBarrelData('R', R+t_fluid, 't', t_fluid, 'parm', parm, 'L', L+2*t_fluid, 'd_p', 3, 'Xi', Xi);
    if strcmp(task.misc.method,'PML')
        pmlLayer = getBarrelData('R', R+2*t_fluid, 't', t_pml, 'parm', parm, 'L', L+4*t_fluid, 'd_p', 3, 'Xi', Xi,'pmlFill', task.msh.pmlFill);
        if parm == 1
            pmlLayer{1}.isPML = [0,0,1];
            pmlLayer{2}.isPML = [0,1,1];
            pmlLayer{3}.isPML = [0,0,1];
            pmlLayer{4}.isPML = [0,0,1];
            pmlLayer{5}.isPML = [0,1,1];
            if task.msh.pmlFill
                pmlLayer{2}.isPML = [0,0,1];
                pmlLayer{5}.isPML = [0,0,1];
            end
        else
            error('Not properly implemented: Needs to make an anulus around special parametrization in the pmlLayer')
            for i = 1:5
                pmlLayer{i}.isPML = [0,0,1];
            end
            for i = 6:9
                pmlLayer{i}.isPML = [0,1,1];
            end
            for i = 10:13
                pmlLayer{i}.isPML = [0,0,1];
            end
            for i = 14:18
                pmlLayer{i}.isPML = [0,1,0];
            end
            for i = 19:22
                pmlLayer{i}.isPML = [0,1,1];
            end
            if task.msh.pmlFill
                for i = 6:9
                    pmlLayer{i}.isPML = [0,0,1];
                end
                for i = 19:22
                    pmlLayer{i}.isPML = [0,1,0];
                end
            end
        end
            
        pmlLayer = insertKnotsInNURBS(pmlLayer,{{[] R/(R+t_fluid) []}, {}, {[] [t_fluid/(L+2*t_fluid), (L+t_fluid)/(L+2*t_fluid)] []}, {[] t_fluid/(R+t_fluid) []}, {}});
        fluid = explodeNURBS(fluid);
        varCol{1}.nurbs = fluid;
        varCol = copySet(varCol,1, 'outer', 'Gamma_a');
        fluid = explodeNURBS([fluid,pmlLayer]);
    end
    fluid = makeUniformNURBSDegree(fluid,degree);
    if parm == 1
        Imap{1} = [R*xiAngle, (R+t_fluid)*xiAngle, (R+t_fluid+t_pml)*xiAngle];
        [fluid,newKnotsIns] = refineNURBSevenly(fluid,(2^(M-1)-1)/R,Imap,0,dirs,true,true);
    else
        Imap{1} = [R*xiAngle, R*1.414057108745095];
        [fluid,newKnotsIns] = refineNURBSevenly(fluid,(2^(M-1)-1)/(R*xiAngle),Imap);
    end
    task.iem.N = min(numel(newKnotsIns{1}{3})+degree(1),9);
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
    varCol{1}.c_xy = c_xy;
end
varCol{1}.chimin = chimin;
varCol{1}.chimax = chimax;
varCol{1}.L_gamma = L;
task.iem.Upsilon = Upsilon;
task.varCol = varCol;