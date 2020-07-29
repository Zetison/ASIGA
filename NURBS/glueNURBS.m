function nurbs = glueNURBS(nurbsCol,dir)
nurbsCol = makeUniformNURBSDegree(nurbsCol);
nurbsCol = makeUniformNURBSDimension(nurbsCol);
noPatches = numel(nurbsCol);
n = zeros(noPatches,1);
d = nurbsCol{1}.d;
d_p = nurbsCol{1}.d_p;
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
        m = size(nurbsCol{1}.coeffs,3);
        coeffs = zeros(d+1,sum(n)-noPatches+1,m);
        p = nurbsCol{1}.degree(1);
        Xi = zeros(1,sum(n)+p+1-noPatches+1);
        if d_p > 1
            Eta = nurbsCol{1}.knots{2};
        end
        Xi(end-p:end) = 1;
    case 'eta'
        m = nurbsCol{1}.number(1);
        coeffs = zeros(d+1,m,sum(n)-noPatches+1);
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
if d_p == 1
    nurbs = createNURBSobject(coeffs,Xi);
else
    nurbs = createNURBSobject(coeffs,{Xi,Eta});
end