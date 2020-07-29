function [nDofs, rednDofs, h_arr, error] = getDataFromFile(filename)


fid = fopen(filename,'r','b');

row = fscanf(fid,'%s %s %s %s\n',4);

formatSpec = '%d %d %f %f';
sizeA = [4 Inf];

A = fscanf(fid,formatSpec,sizeA);

nDofs = A(1,:);
rednDofs = A(2,:);
h_arr = A(3,:);
error = A(4,:);
fclose(fid);