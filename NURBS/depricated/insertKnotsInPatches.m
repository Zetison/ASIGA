function nurbs = insertKnotsInPatches(nurbs,noNewXiKnots,noNewEtaKnots,noNewZetaKnots)

error('Depricated: use insertKnotsInNURBS() instead')
% if ~iscell(nurbs)
%     nurbs = {nurbs};
% end
% for i = 1:numel(nurbs)  
%     switch nurbs{i}.type
%         case '3Dvolume'
%             if numel(noNewXiKnots) > 1
%                 nurbs{i} = insertKnotsInNURBS(nurbs{i},{insertUniform2(nurbs{i}.knots{1}, noNewXiKnots(1)) ...
%                                                         insertUniform2(nurbs{i}.knots{2}, noNewXiKnots(2)) ...
%                                                         insertUniform2(nurbs{i}.knots{3}, noNewXiKnots(3))}); 
%             else
%                 nurbs{i} = insertKnotsInNURBS(nurbs{i},{insertUniform2(nurbs{i}.knots{1}, noNewXiKnots) ...
%                                                         insertUniform2(nurbs{i}.knots{2}, noNewEtaKnots) ...
%                                                         insertUniform2(nurbs{i}.knots{3}, noNewZetaKnots)}); 
%             end
%         case {'3Dsurface','2Dsurface'}
%             if numel(noNewXiKnots) > 1
%                 nurbs{i} = insertKnotsInNURBS(nurbs{i},{insertUniform2(nurbs{i}.knots{1}, noNewXiKnots(1)) ...
%                                                         insertUniform2(nurbs{i}.knots{2}, noNewXiKnots(2))}); 
%             else
%                 nurbs{i} = insertKnotsInNURBS(nurbs{i},{insertUniform2(nurbs{i}.knots{1}, noNewXiKnots) ...
%                                                         insertUniform2(nurbs{i}.knots{2}, noNewEtaKnots)}); 
%             end
%     end
% end