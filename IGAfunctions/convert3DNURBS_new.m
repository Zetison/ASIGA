function varCol = convert3DNURBS_new(nurbs, varCol)

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

        varCol.weights = weights;
        varCol.controlPts = controlPts;


        varCol.nurbs = nurbs;
        varCol.noCtrlPts = nurbs.number(1)*nurbs.number(2)*nurbs.number(3);
        varCol.noDofs = varCol.dimension*varCol.noCtrlPts;
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

        varCol.weights = weights;
        varCol.controlPts = controlPts;


        varCol.nurbs = nurbs;
        varCol.noCtrlPts = nurbs.number(1)*nurbs.number(2);
        varCol.noDofs = varCol.dimension*varCol.noCtrlPts;
    case '1Dnurbs'
        p = nurbs.degree;

        Xi   = nurbs.knots{1};

        n = length(Xi)-p-1;

        weights = reshape(nurbs.coeffs(2,:),n,1);

        controlPts = nurbs.coeffs(1,:);
        
        varCol.weights = weights;
        varCol.controlPts = controlPts;


        varCol.nurbs = nurbs;
        varCol.noCtrlPts = nurbs.number;
        varCol.noDofs = varCol.dimension*varCol.noCtrlPts;
end