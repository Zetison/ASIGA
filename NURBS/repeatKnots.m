function varCol = repeatKnots(varCol,coreMethod)
nurbsPatches = varCol.nurbs;
if ~iscell(nurbsPatches)
    nurbsPatches = {nurbsPatches};
end
for patch = 1:numel(nurbsPatches)
    nurbs = nurbsPatches{patch};
    switch coreMethod
        case {'linear_FEM', 'h_FEM', 'hp_FEM', 'C0_IGA', 'SEM'}
            if strcmp(nurbs.type, '3Dvolume')
                p_xi = nurbs.degree(1);
                p_eta = nurbs.degree(2);
                p_zeta = nurbs.degree(3);

                Xi = nurbs.knots{1};
                Eta = nurbs.knots{2};
                Zeta = nurbs.knots{3};

                for i = 1:length(Xi)
                    mm = length(find(Xi == Xi(i)));    
                    nurbs = insertKnotsInNURBS(nurbs,{Xi(i)*ones(p_xi-mm,1) [] []});
                end

                for j = 1:length(Eta)
                    mm = length(find(Eta == Eta(j)));    
                    nurbs = insertKnotsInNURBS(nurbs,{[] Eta(j)*ones(p_eta-mm,1) []});
                end

                for k = 1:length(Zeta)
                    mm = length(find(Zeta == Zeta(k)));    
                    nurbs = insertKnotsInNURBS(nurbs,{[] [] Zeta(k)*ones(p_zeta-mm,1)});
                end
            elseif strcmp(nurbs.type, '3Dsurface')
                p_xi = nurbs.degree(1);
                p_eta = nurbs.degree(2);

                Xi = nurbs.knots{1};
                Eta = nurbs.knots{2};


                for i = 1:length(Xi)
                    mm = length(find(Xi == Xi(i)));    
                    nurbs = insertKnotsInNURBS(nurbs,{Xi(i)*ones(p_xi-mm,1) []});
                end

                for j = 1:length(Eta)
                    mm = length(find(Eta == Eta(j)));    
                    nurbs = insertKnotsInNURBS(nurbs,{[] Eta(j)*ones(p_eta-mm,1)});
                end



            end
    end
    nurbsPatches{patch} = nurbs;
end
varCol.nurbs = nurbsPatches;