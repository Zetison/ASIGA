function nurbs = glueNURBS(nurbsCol,dir)
nurbsCol = makeUniformNURBSDegree(nurbsCol);
nurbsCol = makeUniformNURBSDimension(nurbsCol);
noPatches = numel(nurbsCol);
n = zeros(noPatches,1);
for i = 1:noPatches
    n(i) = nurbsCol{i}.number(dir);
end

dimension = size(nurbsCol{1}.coeffs);
dimension(dir+1) = sum(n)-noPatches+1;
coeffs = zeros(dimension);
p = nurbsCol{1}.degree(dir);
Xi = zeros(1,sum(n)+p+1-noPatches+1);
Xi(end-p:end) = 1;

counter = 1;
counter2 = 2;
for i = 1:noPatches
    coeffs = subasgnArr(coeffs,nurbsCol{i}.coeffs,counter:counter+n(i)-1,dir+1);
    Xi(counter2:counter2-1+n(i)-1) = nurbsCol{i}.knots{dir}(2:end-p-1)/noPatches + (i-1)/noPatches;    
    counter  = counter  + n(i) - 1;
    counter2 = counter2 + n(i) - 1;
end
knots = nurbsCol{1}.knots;
knots{dir} = Xi;
nurbs = createNURBSobject(coeffs,knots);