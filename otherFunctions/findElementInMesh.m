function idx = findElementInMesh(xi, elRangeXi)

if xi == 0
    idx = 1;
    return
end
i = size(elRangeXi,1);
while xi < elRangeXi(i,1)
    i = i - 1;
end
idx = i;