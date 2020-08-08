function refinednurbs = refineNURBSevenly(nurbs,n,Imap,noExtraEvalPts)
if nargin < 3
    Imap = {};
end
if nargin < 4
    noExtraEvalPts = 1;
end
for i = 1:numel(Imap)
    Imap{i} = fliplr(sort(Imap{i}));
end
Eps = 1e-10;
refinednurbs = nurbs;
I_maxRecord = [];
for patch = 1:numel(nurbs)
    d_p = nurbs{patch}.d_p;
    if d_p < 2 || d_p > 3
        error('Not implemented for this case')
    end
    for dir = 1:d_p
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
                I = NURBSarcLength(nurbs{patch},uniqueXi(j),uniqueXi(j+1),parm_pts(i,:),dir);
                if I_max < I
                    I_max = I;
                end
            end
            for i = 1:numel(Imap)
                if any(abs(Imap{i}-I_max) < Eps)
                    I_max = Imap{i}(1);
                    break
                end
            end
            if ~any(abs(I_maxRecord-I_max) < Eps)
                I_maxRecord = [I_maxRecord, I_max];
            end
            newXiKnots = linspace2(uniqueXi(j),uniqueXi(j+1),round(I_max*n));
            newKnots = cell(1,d_p);
            newKnots{dir} = newXiKnots;
            refinednurbs(patch) = insertKnotsInNURBS(refinednurbs(patch),newKnots);
        end
    end
end
% I_maxRecord