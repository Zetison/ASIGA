function halfSizeOfDataSet(oldFile, newFile)


fid = fopen(oldFile,'r','b');

row = fscanf(fid,'%s %s\n',2);

formatSpec = '%f %f';
sizeA = [2 Inf];

A = fscanf(fid,formatSpec,sizeA);

x = A(1,:);
y = A(2,:);

x = x(1:2:end);
y = y(1:2:end);

fclose(fid);

fid = fopen(newFile,'wt+','b');
fprintf(fid,'k\t\t\tF_k\n');
for j = 1:length(x)
    fprintf(fid,'%1.15f\t%1.15f\n',x(j),abs(y(j)));
end
fclose(fid);


