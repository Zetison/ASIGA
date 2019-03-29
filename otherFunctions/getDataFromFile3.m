function [x, y] = getDataFromFile3(filename)


fid = fopen(filename,'r','b');
if fid == -1
    x = NaN;
    y = NaN;
    return
end

row = fscanf(fid,'%s %s\n',2);

formatSpec = '%f %f';
sizeA = [2 Inf];

A = fscanf(fid,formatSpec,sizeA);

x = A(1,:);
y = A(2,:);
fclose(fid);