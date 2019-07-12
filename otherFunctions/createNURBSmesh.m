function [varCol, fluid, solid, fluid_i] = createNURBSmesh(varCol, parms, model, M, degree)

solid = NaN;
fluid_i = NaN;
switch model
    case 'EL'
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_Ellipsoid(varCol,parms, M, degree, model); 
    case {'M1','M1_P'}
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_Model1(varCol,parms, M, degree);
        error('Clean up this subRoutine')
    case {'M3','M3_P'}
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_Model3(varCol,parms, M, degree);
    case {'M3','M3_P'}
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_Model2(varCol,parms, M, degree);
    case {'PH','PH_P'}
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_PH(varCol,parms, M, degree);
    case {'M4','M4_P'}
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_Model4(varCol,parms, M, degree);
    case {'Torus'}
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_Torus(varCol,parms, M, degree);
    case {'Cube','Cube_P'}
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_Cube(varCol,parms, M, degree);
    case {'M5A','M5B','M5A_P','M5B_P'}
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_Model5(varCol,parms, M, degree);
        error('Clean up this subRoutine')
    case {'MS','MS_P'}
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_MockShell(varCol,parms, M, degree);
    case 'TAP'
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_TAP(varCol,parms, M, degree);
    case {'BA','BA_P'}
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_Barrel(varCol,parms, M, degree);
        error('Clean up this subRoutine')
    case {'BC','BC_P'}
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_BC(varCol,parms, M, degree);
    case {'BCA','BCA_P'}
        [varCol, fluid, solid, fluid_i] = createNURBSmesh_BCA(varCol,parms, M, degree);
    otherwise
        if varCol.isSphericalShell
            [varCol, fluid, solid, fluid_i] = createNURBSmesh_Ellipsoid(varCol, parms, M, degree, model);
        end
end


coreMethod = varCol.coreMethod;
useSolidDomain = varCol.useSolidDomain;
useInnerFluidDomain = varCol.useInnerFluidDomain;

fluid = repeatKnots(fluid,coreMethod);
if useSolidDomain
    solid = repeatKnots(solid,coreMethod);
end
if useInnerFluidDomain
    fluid_i = repeatKnots(fluid_i,coreMethod);
end
fluid = degenerateIGAtoFEM_Grev(fluid,coreMethod);
if useSolidDomain
    solid = degenerateIGAtoFEM_Grev(solid,coreMethod);
end
if useInnerFluidDomain
    fluid_i = degenerateIGAtoFEM_Grev(fluid_i,coreMethod);
end
