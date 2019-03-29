function varCol = convertNURBS(nurbsPatches, varCol)

if ~iscell(nurbsPatches)
    nurbsPatches = {nurbsPatches};
end
noPatches = numel(nurbsPatches);
patches = cell(1,noPatches);
for patch = 1:numel(nurbsPatches)
    nurbs = nurbsPatches{patch};
    switch nurbs.type
        case '3Dvolume'
            p = nurbs.degree(1);
            q = nurbs.degree(2);
            r = nurbs.degree(3);

            Xi   = nurbs.knots{1};
            Eta  = nurbs.knots{2};
            Zeta = nurbs.knots{3};

            n = length(Xi)-p-1;
            m = length(Eta)-q-1;
            l = length(Zeta)-r-1;

            weights = reshape(nurbs.coeffs(4,:,:),n*m*l,1);

            controlPts = zeros(n*m*l,3);

            count = 0;
            for k=1:l
                for j=1:m
                    controlPts(n*count+1:n*(count+1),:) = nurbs.coeffs(1:3,:,j,k)';
                    count = count+1;
                end
            end

            patches{patch}.weights = weights;
            patches{patch}.controlPts = controlPts;


            patches{patch}.nurbs = nurbs;
            patches{patch}.noCtrlPts = nurbs.number(1)*nurbs.number(2)*nurbs.number(3);
            patches{patch}.noDofs = varCol.dimension*patches{patch}.noCtrlPts;
        case '3Dsurface'
            p = nurbs.degree(1);
            q = nurbs.degree(2);

            Xi   = nurbs.knots{1};
            Eta  = nurbs.knots{2};

            n = length(Xi)-p-1;
            m = length(Eta)-q-1;

            weights = reshape(nurbs.coeffs(4,:,:),n*m,1);

            controlPts = zeros(n*m,3);

            count = 0;
            for j=1:m
                controlPts(n*count+1:n*(count+1),:) = nurbs.coeffs(1:3,:,j)';
                count = count+1;
            end

            patches{patch}.weights = weights;
            patches{patch}.controlPts = controlPts;


            patches{patch}.nurbs = nurbs;
            patches{patch}.noCtrlPts = nurbs.number(1)*nurbs.number(2);
            patches{patch}.noDofs = varCol.dimension*patches{patch}.noCtrlPts;
        case '2Dsurface'
            p = nurbs.degree(1);
            q = nurbs.degree(2);

            Xi   = nurbs.knots{1};
            Eta  = nurbs.knots{2};

            n = length(Xi)-p-1;
            m = length(Eta)-q-1;

            weights = reshape(nurbs.coeffs(3,:,:),n*m,1);

            controlPts = zeros(n*m,2);

            count = 0;
            for j=1:m
                controlPts(n*count+1:n*(count+1),:) = nurbs.coeffs(1:2,:,j)';
                count = count+1;
            end

            patches{patch}.weights = weights;
            patches{patch}.controlPts = controlPts;


            patches{patch}.nurbs = nurbs;
            patches{patch}.noCtrlPts = nurbs.number(1)*nurbs.number(2);
            patches{patch}.noDofs = varCol.dimension*patches{patch}.noCtrlPts;
        case '2Dcurve'
            p = nurbs.degree;

            Xi = nurbs.knots;

            n = length(Xi)-p-1;

            weights = reshape(nurbs.coeffs(3,:),n,1);

            controlPts = nurbs.coeffs(1:2,:)';

            patches{patch}.weights = weights;
            patches{patch}.controlPts = controlPts;


            patches{patch}.nurbs = nurbs;
            patches{patch}.noCtrlPts = nurbs.number;
            patches{patch}.noDofs = varCol.dimension*patches{patch}.noCtrlPts;
        case '1Dnurbs'
            p = nurbs.degree;

            Xi = nurbs.knots;

            n = length(Xi)-p-1;

            weights = reshape(nurbs.coeffs(2,:),n,1);

            controlPts = nurbs.coeffs(1,:)';

            patches{patch}.weights = weights;
            patches{patch}.controlPts = controlPts;


            patches{patch}.nurbs = nurbs;
            patches{patch}.noCtrlPts = nurbs.number;
            patches{patch}.noDofs = varCol.dimension*patches{patch}.noCtrlPts;
    end
end
varCol.patches = patches;