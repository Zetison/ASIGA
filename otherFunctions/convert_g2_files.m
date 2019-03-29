function convert_g2_files(inputFileName,outputFileName)


fid = fopen(inputFileName,'r','b');

A = fscanf(fid,'%d\n',4);
A = fscanf(fid,'%d\n',2);


A = fscanf(fid,'%d\n',2);
n = A(1);
p = A(2)-1;
Xi = fscanf(fid,'%f\n',n+p+1)';

A = fscanf(fid,'%d\n',2);
m = A(1);
q = A(2)-1;
Eta = fscanf(fid,'%f\n',m+q+1)';

A = fscanf(fid,'%d\n',2);
l = A(1);
r = A(2)-1;
Zeta = fscanf(fid,'%f\n',l+r+1)';

controlPts = zeros(4,n,m,l);
for k = 1:l
    for j = 1:m
        for i = 1:n
            controlPts(:,i,j,k) = fscanf(fid,'%f\n',4)';
        end
    end
end
            

for i = 1:n
    for j = 1:m
        for k = 1:l
            controlPts(1:3,i,j,k) = controlPts(1:3,i,j,k)/controlPts(4,i,j,k);
        end
    end
end

fclose(fid);

fid = fopen(outputFileName,'wt+','b');
fprintf(fid,'Xi = [');
for i = 1:length(Xi)
    fprintf(fid,'%f ',Xi(i));
end
fprintf(fid,'];\n');

fprintf(fid,'Eta = [');
for i = 1:length(Eta)
    fprintf(fid,'%f ',Eta(i));
end
fprintf(fid,'];\n');

fprintf(fid,'Zeta = [');
for i = 1:length(Zeta)
    fprintf(fid,'%f ',Zeta(i));
end
fprintf(fid,'];\n');

fprintf(fid,'\n');

fprintf(fid,'controlPts = zeros(4,%d,%d,%d);\n',n,m,l);

fprintf(fid,'\n');
for k = 1:l
    for j = 1:m
        for i = 1:n
            fprintf(fid,'controlPts(:,%d,%d,%d) = [  \t%1.15f      \t%1.15f       \t%1.15f    \t%1.15f           \t];\n',i,j,k,controlPts(:,i,j,k));
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'solid = createNURBSobject(controlPts,{Xi, Eta, Zeta});');


fclose(fid);