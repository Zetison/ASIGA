function nurbs = insertKnotsInPatches(nurbs,noNewXiKnots,noNewEtaKnots,noNewZetaKnots)


if ~iscell(nurbs)
    nurbs = {nurbs};
end
for i = 1:numel(nurbs)  
    switch nurbs{i}.type
        case '3Dvolume'
            nurbs{i} = insertKnotsInNURBS(nurbs{i},{insertUniform2(nurbs{i}.knots{1}, noNewXiKnots) ...
                                                    insertUniform2(nurbs{i}.knots{2}, noNewEtaKnots) ...
                                                    insertUniform2(nurbs{i}.knots{3}, noNewZetaKnots)}); 
        case {'3Dsurface','2Dsurface'}
            nurbs{i} = insertKnotsInNURBS(nurbs{i},{insertUniform2(nurbs{i}.knots{1}, noNewXiKnots) ...
                                                    insertUniform2(nurbs{i}.knots{2}, noNewEtaKnots)}); 
    end
end