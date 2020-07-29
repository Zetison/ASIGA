function varCol = convertNURBS(varCol)
if ~isfield(varCol,'nurbs')
    return
end
nurbsPatches = varCol.nurbs;
if ~iscell(nurbsPatches)
    nurbsPatches = {nurbsPatches};
end
noPatches = numel(nurbsPatches);
patches = cell(1,noPatches);
for patch = 1:numel(nurbsPatches)
    nurbs = nurbsPatches{patch};
    d = nurbs.d;
    switch nurbs.d_p
        case 3
            p = nurbs.degree(1);
            q = nurbs.degree(2);
            r = nurbs.degree(3);

            Xi   = nurbs.knots{1};
            Eta  = nurbs.knots{2};
            Zeta = nurbs.knots{3};

            n = length(Xi)-p-1;
            m = length(Eta)-q-1;
            l = length(Zeta)-r-1;

            weights = reshape(nurbs.coeffs(d+1,:,:),n*m*l,1);

            controlPts = zeros(n*m*l,d);

            count = 0;
            for k=1:l
                for j=1:m
                    controlPts(n*count+1:n*(count+1),:) = nurbs.coeffs(1:d,:,j,k)';
                    count = count+1;
                end
            end

            patches{patch}.weights = weights;
            patches{patch}.controlPts = controlPts;

            patches{patch}.noCtrlPts = nurbs.number(1)*nurbs.number(2)*nurbs.number(3);
        case 2
            p = nurbs.degree(1);
            q = nurbs.degree(2);

            Xi   = nurbs.knots{1};
            Eta  = nurbs.knots{2};

            n = length(Xi)-p-1;
            m = length(Eta)-q-1;

            weights = reshape(nurbs.coeffs(d+1,:,:),n*m,1);

            controlPts = zeros(n*m,d);

            count = 0;
            for j=1:m
                controlPts(n*count+1:n*(count+1),:) = nurbs.coeffs(1:d,:,j)';
                count = count+1;
            end

            patches{patch}.weights = weights;
            patches{patch}.controlPts = controlPts;

            patches{patch}.noCtrlPts = nurbs.number(1)*nurbs.number(2);
        case 1
            p = nurbs.degree;

            Xi = nurbs.knots;

            n = length(Xi)-p-1;

            weights = reshape(nurbs.coeffs(d+1,:),n,1);

            controlPts = nurbs.coeffs(1:d,:)';

            patches{patch}.weights = weights;
            patches{patch}.controlPts = controlPts;

            patches{patch}.noCtrlPts = nurbs.number;
    end
    patches{patch}.nurbs = nurbs;
    patches{patch}.noDofs = varCol.dimension*patches{patch}.noCtrlPts;
end
varCol.patches = patches;