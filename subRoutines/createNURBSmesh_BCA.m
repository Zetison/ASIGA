function [varCol, fluid, solid, fluid_i] = createNURBSmesh_BCA(varCol, parms, M, degree)

solid = NaN;
fluid_i = NaN;


names = fieldnames(parms);
for j = 1:numel(names)
    name = names{j};
    eval([name, ' = parms.(names{j});']);
end
alpha = parms.alpha; % resolve matlab bug
beta = parms.beta; % resolve matlab bug

if varCol.boundaryMethod
    noNewXiKnots = 2^(M-1)-1;
    noNewEtaKnots = noNewXiKnots;
    if false
        nurbsColApprox = getBCAData(a, b, L, g2, g3, alpha, beta, c, s, degree);
    else
        load(['NURBSgeometries/BCAdata/BeTSSi_BCA_p' num2str(degree)])
        fluid = nurbs;
    end
    
    fluid = insertKnotsInPatches(fluid,noNewXiKnots,noNewEtaKnots);
    varCol.patchTop = getPatchTopology(fluid);
end
L_gamma = L+a+g2+g3;

varCol.L_gamma = L_gamma;
