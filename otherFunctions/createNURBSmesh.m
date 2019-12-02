function varCol = createNURBSmesh(varCol, parms, model, M, degree)

switch model
    case 'EL'
        varCol = createNURBSmesh_Ellipsoid(varCol,parms, M, degree, model); 
    case {'M1','M1_P'}
        varCol = createNURBSmesh_Model1(varCol,parms, M, degree);
    case {'M2','M2_P'}
        varCol = createNURBSmesh_Model2(varCol,parms, M, degree);
    case {'M3','M3_P'}
        varCol = createNURBSmesh_Model3(varCol,parms, M, degree);
    case 'M31'
        varCol = createNURBSmesh_Model3(varCol,parms, M, degree);
    case 'S15'
        varCol = createNURBSmesh_S15(varCol,parms, M, degree);
    case {'PH','PH_P'}
        varCol = createNURBSmesh_PH(varCol,parms, M, degree);
    case {'M4','M4_P'}
        varCol = createNURBSmesh_Model4(varCol,parms, M, degree);
    case {'Torus'}
        varCol = createNURBSmesh_Torus(varCol,parms, M, degree);
    case {'Cube','Cube_P'}
        varCol = createNURBSmesh_Cube(varCol,parms, M, degree);
    case {'M5A','M5B','M5A_P','M5B_P'}
        varCol = createNURBSmesh_Model5(varCol,parms, M, degree);
    case {'MS','MS_P'}
        varCol = createNURBSmesh_MockShell(varCol,parms, M, degree);
    case {'Shirron'}
        varCol = createNURBSmesh_Shirron(varCol,parms, M, degree);
    case 'TAP'
        varCol = createNURBSmesh_TAP(varCol,parms, M, degree);
    case {'Barrel','Barrel_P'}
        varCol = createNURBSmesh_Barrel(varCol,parms, M, degree);
    case {'BC','BC_P'}
        varCol = createNURBSmesh_BC(varCol,parms, M, degree);
    case {'BCA','BCA_P'}
        varCol = createNURBSmesh_BCA(varCol,parms, M, degree);
    otherwise
        if varCol{1}.isSphericalShell
            varCol = createNURBSmesh_Ellipsoid(varCol, parms, M, degree, model);
        end
end

for i = 1:numel(varCol) % assume coreMethod to be the same in all domains
    varCol{i} = repeatKnots(varCol{i},varCol{1}.coreMethod);
    varCol{i} = degenerateIGAtoFEM_Grev(varCol{i},varCol{1}.coreMethod);
end
