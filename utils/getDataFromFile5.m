function [x, y, z] = getDataFromFile5(filename)


fid = fopen(filename,'r','b');

row = fscanf(fid,'%s\n',4);

formatSpec = '%f %f %f %f';
sizeA = [4 Inf];

A = fscanf(fid,formatSpec,sizeA);

x = A(1,:);
y = A(2,:);
z = A(3,:);
fclose(fid);