function [nurbs, maxLengths] = autoRefineNURBS(nurbs,connection,h_max,dirs)

d_p = nurbs{1}.d_p;
if nargin < 4
    dirs = 1:d_p;
end



noPatches = numel(nurbs);
patchChecked = false(noPatches,d_p);
maxLengths = zeros(noPatches,d_p);
topologyMap = cell(1,noPatches);
for patch = 1:noPatches
    topologyMap{patch}(2*d_p) = struct();
    for i = 1:numel(connection)
        if connection{i}.Attributes.master == patch
            midx = connection{i}.Attributes.midx;
            slave = connection{i}.Attributes.slave;
            sidx = connection{i}.Attributes.sidx;
            orient = connection{i}.Attributes.orient;
            topologyMap{patch}(midx).slave = slave;
            topologyMap{patch}(midx).sidx = sidx;
            topologyMap{patch}(midx).orient = orient;
        end
    end
end
% I = edgeLengths(nurbs,dirs);
maxLengths = zeros(noPatches,d_p);
% for patch = 1:noPatches
%     for i = 1:d_p
%         maxLengths(patch,i) = max(I(i,:));
%     end
% end


for patch = 1:noPatches
    [maxLengths, patchChecked] = searchNURBS(nurbs,topologyMap,maxLengths, patchChecked,dirs,patch);
end
for patch = 1:noPatches
    for i = 1:numel(dirs)
        nurbs = insertKnotsInNURBS(nurbs,round(maxLengths(patch,i)/h_max)-1);
    end
end

function [maxLengths, patchChecked] = searchNURBS(nurbs,topologyMap,maxLengths, patchChecked,dirs,patch)

if numel(dirs) == 1
    dir = topologyMap{patch}(midx).slave;
    [maxLengths, patchChecked] = searchNURBS(nurbs,topologyMap,maxLengths,patchChecked,dir,patch);
else
    for dir = dirs
        [maxLengths, patchChecked] = searchNURBS(nurbs,topologyMap,maxLengths,patchChecked,dir,patch);
    end
end
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













% function [maxArcLength, patchRefined] = calcArcLengths(nurbs,patch,dir,maxArcLength, patchRefined,type,refLength)
% Eps = 1e-10;
% I_maxRecord = [];
% maxArcLength = cell(1,noPatches);
% d_p = nurbs{patch}.d_p;
% newKnots = cell(1,d_p);
% maxArcLength{patch} = cell(1,d_p);
% if d_p == 1 || d_p >= 4
%     error('Not implemented for this case')
% end
% uniqueXi = unique(nurbs{patch}.knots{dir});
% dir2 = setdiff(1:d_p,dir);
% 
% knots = nurbs{patch}.knots(dir2);
% for j = 1:numel(dir2)
%     knots{j} = insertUniform(unique(knots{j}), noExtraEvalPts);
% end
% if d_p == 2
%     uniqueEta = knots{1};
% elseif d_p == 3
%     uniqueEta = [kron(knots{1},ones(numel(knots{2}),1)), kron(ones(numel(knots{1}),1),knots{2})];
% end
% parm_pts = NaN(size(uniqueEta,1),d_p);
% parm_pts(:,dir2) = uniqueEta;
% for j = 1:numel(uniqueXi)-1
%     I_max = -Inf;
%     for i = 1:size(parm_pts,1)
%         I = NURBSarcLength(nurbs{patch},uniqueXi(j),uniqueXi(j+1),parm_pts(i,:),dir,true);
%         if I_max < I
%             I_max = I;
%         end
%     end
%     maxArcLength(j) = I_max;
% end
% maxArcLength = max(maxArcLength);