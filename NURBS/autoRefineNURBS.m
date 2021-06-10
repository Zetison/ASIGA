function [nurbs, newKnotsIns] = autoRefineNURBS(nurbs,refLength,n,dirs)

if nargin < 4
    d_p = nurbs{1}.d_p;
    dirs = 1:d_p;
end


Eps = 1e-10;
I_maxRecord = [];

noPatches = numel(nurbs);
patchRefined = false(noPatches,3);
newKnotsIns = cell(1,noPatches);
for patch = 1:noPatches
    [newKnotsIns, patchRefined] = searchNURBS(nurbs,newKnotsIns, patchRefined,dirs,'findLength',refLength,n);
end
nurbs = insertKnotsInNURBS(nurbs,newKnotsIns);

function [maxArcLength, patchRefined] = calcArcLengths(nurbs,maxArcLength, patchRefined,dirs,type,refLength,n)
maxArcLength = cell(1,noPatches);
for patch = 1:noPatches
    d_p = nurbs{patch}.d_p;
    newKnots = cell(1,d_p);
    maxArcLength{patch} = cell(1,d_p);
    if d_p == 1 || d_p >= 4
        error('Not implemented for this case')
    end
    for dir = dirs
        uniqueXi = unique(nurbs{patch}.knots{dir});
        dir2 = setdiff(1:d_p,dir);

        knots = nurbs{patch}.knots(dir2);
        for j = 1:numel(dir2)
            knots{j} = insertUniform(unique(knots{j}), noExtraEvalPts);
        end
        if d_p == 2
            uniqueEta = knots{1};
        elseif d_p == 3
            uniqueEta = [kron(knots{1},ones(numel(knots{2}),1)), kron(ones(numel(knots{1}),1),knots{2})];
        end
        parm_pts = NaN(size(uniqueEta,1),d_p);
        parm_pts(:,dir2) = uniqueEta;
        for j = 1:numel(uniqueXi)-1
            I_max = -Inf;
            for i = 1:size(parm_pts,1)
                I = NURBSarcLength(nurbs{patch},uniqueXi(j),uniqueXi(j+1),parm_pts(i,:),dir,true);
                if I_max < I
                    I_max = I;
                end
            end
            maxArcLength{patch}{dir}(j) = I_max;
        end
    end
end
maxArcLength{patch} = newKnots;

function [newKnotsIns, patchRefined] = searchNURBS(nurbs,newKnotsIns, patchRefined,dirs,type,refLength,n)
% if nargin < 3
%     Imap = {};
% end
% 
% if nargin < 4
%     noExtraEvalPts = 1;
% end
% if nargin < 5
%     d_p = nurbs{1}.d_p;
%     dirs = 1:d_p;
% end
% if nargin < 6
%     sortImap = true;
% end
% if sortImap
%     for i = 1:numel(Imap)
%         Imap{i} = fliplr(sort(Imap{i}));
%     end
% end
% for patch = 1:noPatches
%     d_p = nurbs{patch}.d_p;
%     if d_p < 2 || d_p > 3
%         error('Not implemented for this case')
%     end
%     newKnots = cell(1,d_p);
%     for dir = dirs
%         uniqueXi = unique(nurbs{patch}.knots{dir});
%         dir2 = setdiff(1:d_p,dir);
%         
%         knots = nurbs{patch}.knots(dir2);
%         for j = 1:numel(dir2)
%             knots{j} = insertUniform(unique(knots{j}), noExtraEvalPts);
%         end
%         if d_p == 2
%             uniqueEta = knots{1};
%         elseif d_p == 3
%             uniqueEta = [kron(knots{1},ones(numel(knots{2}),1)), kron(ones(numel(knots{1}),1),knots{2})];
%         end
%         parm_pts = NaN(size(uniqueEta,1),d_p);
%         parm_pts(:,dir2) = uniqueEta;
%         newXiKnots = [];
%         for j = 1:numel(uniqueXi)-1
%             I_max = -Inf;
%             for i = 1:size(parm_pts,1)
%                 I = NURBSarcLength(nurbs{patch},uniqueXi(j),uniqueXi(j+1),parm_pts(i,:),dir,true);
%                 if I_max < I
%                     I_max = I;
%                 end
%             end
%             for i = 1:numel(Imap)
%                 if any(abs(Imap{i}-I_max) < Eps)
%                     I_max = Imap{i}(1);
%                     break
%                 end
%             end
%             if ~any(abs(I_maxRecord-I_max) < Eps)
%                 I_maxRecord = [I_maxRecord, I_max];
%             end
%             newXiKnots = [newXiKnots, linspace2(uniqueXi(j),uniqueXi(j+1),round(I_max*n))];
%         end
%         newKnots{dir} = newXiKnots;
%     end
%     newKnotsIns{patch} = newKnots;
% end
% nurbs = insertKnotsInNURBS(nurbs,newKnotsIns);
% % I_maxRecord