function [x, y, z, w] = getDataFromFile2(filename)


fid = fopen(filename,'r','b');

row = fscanf(fid,'%s %s %s %s\n',4);

formatSpec = '%f %f %f %f';
sizeA = [4 Inf];

A = fscanf(fid,formatSpec,sizeA);

x = A(1,:);
y = A(2,:);
z = A(3,:);
w = A(4,:);
fclose(fid);