function varCol = createNURBSmesh_BCA(varCol, M, degree)


if varCol{1}.boundaryMethod
    noNewXiKnots = 2^(M-1)-1;
    noNewEtaKnots = noNewXiKnots;
    if false
        nurbsColApprox = getBCAData(a, b, L, g2, g3, alpha, beta, c, s, degree);
    else
        load(['NURBSgeometries/BCAdata/BeTSSi_BCA_p' num2str(degree)])
        fluid = nurbs;
    end
    
    varCol{1}.patchTop = getPatchTopology(fluid);
    fluid = insertKnotsInPatches(fluid,noNewXiKnots,noNewEtaKnots);
end
L_gamma = L+a+g2+g3;

varCol{1}.L_gamma = L_gamma;
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
