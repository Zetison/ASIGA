function A = getDataFromFile_42(filename)


fid = fopen(filename,'r','b');

row = fgets(fid);

formatSpec = '%f %f';
sizeA = [2 Inf];

A = fscanf(fid,formatSpec,sizeA);

x = A(1,:);
y = A(2,:);
fclose(fid);