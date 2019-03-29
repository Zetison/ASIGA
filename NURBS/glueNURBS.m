function nurbs = glueNURBS(nurbsCol,dir)

noPatches = numel(nurbsCol);
n = zeros(noPatches,1);
for i = 1:noPatches
    switch dir
        case 'xi'
            n(i) = nurbsCol{i}.number(1);
        case 'eta'
            n(i) = nurbsCol{i}.number(2);
    end
end
switch dir
    case 'xi'
        m = nurbsCol{1}.number(2);
        coeffs = zeros(4,sum(n)-noPatches+1,m);
        p = nurbsCol{1}.degree(1);
        Xi = zeros(1,sum(n)+p+1-noPatches+1);
        Eta = nurbsCol{1}.knots{2};
        Xi(end-p:end) = 1;
    case 'eta'
        m = nurbsCol{1}.number(1);
        coeffs = zeros(4,m,sum(n)-noPatches+1);
        p = nurbsCol{1}.degree(2);
        Xi = nurbsCol{1}.knots{1};
        Eta = zeros(1,sum(n)+p+1-noPatches+1);
        Eta(end-p:end) = 1;
end
counter = 1;
counter2 = 2;
for i = 1:noPatches
    switch dir
        case 'xi'
            coeffs(:,counter:counter+n(i)-1,:) = nurbsCol{i}.coeffs;
            Xi(counter2:counter2-1+n(i)-1) = nurbsCol{i}.knots{1}(2:end-p-1)/noPatches + (i-1)/noPatches;
        case 'eta'
            coeffs(:,:,counter:counter+n(i)-1) = nurbsCol{i}.coeffs;
            Eta(counter2:counter2-1+n(i)-1) = nurbsCol{i}.knots{2}(2:end-p-1)/noPatches + (i-1)/noPatches;
    end
    counter  = counter  + n(i) - 1;
    counter2 = counter2 + n(i) - 1;
end
nurbs = createNURBSobject(coeffs,{Xi,Eta});