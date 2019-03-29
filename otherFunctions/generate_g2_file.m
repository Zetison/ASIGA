function generate_g2_file(solid,outputFileName)

n = solid.number(1);
m = solid.number(2);
l = solid.number(3);

p = solid.degree(1);
q = solid.degree(2);
r = solid.degree(3);

Xi = solid.knots{1};
Eta = solid.knots{2};
Zeta = solid.knots{3};

controlPts = solid.coeffs;      

for i = 1:n
    for j = 1:m
        for k = 1:l
            controlPts(1:3,i,j,k) = controlPts(1:3,i,j,k)*controlPts(4,i,j,k);
        end
    end
end

fid = fopen(outputFileName,'wt+','b');
fprintf(fid,'700 1 0 0\n');
fprintf(fid,'3 1\n');

fprintf(fid,'%d %d\n',n,p+1);
for i = 1:length(Xi)
    fprintf(fid, '%g ',Xi(i));
end
fprintf(fid, '\n');

fprintf(fid,'%d %d\n',m,q+1);
for j = 1:length(Eta)
    fprintf(fid, '%g ',Eta(j));
end
fprintf(fid, '\n');

fprintf(fid,'%d %d\n',l,r+1);
for k = 1:length(Zeta)
    fprintf(fid, '%g ',Zeta(k));
end
fprintf(fid, '\n');



for k = 1:l
    for j = 1:m
        for i = 1:n
            fprintf(fid, '%.15g %.15g %.15g %.15g\n', controlPts(:,i,j,k));
        end
    end
end
      

fclose(fid);