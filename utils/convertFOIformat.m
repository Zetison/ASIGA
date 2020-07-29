
f = 5;
fid = fopen(['plotData/BeTSSi_M1_HWBC_BI/M1_HWBC_BI_0_0p' num2str(f) '.txt'],'r','b');

row1 = fscanf(fid,'%s\n',1);
row23 = fscanf(fid,'%s\n',1);
row2 = fscanf(fid,'%s\n',1);
row3 = fscanf(fid,'%s\n',2);
row4 = fscanf(fid,'%s\n',1);
row5 = fscanf(fid,'%s\n',1);
row6 = fscanf(fid,'%s\n',1);

formatSpec = '%f, %f, %f';
sizeA = [3 Inf];

A = fscanf(fid,formatSpec,sizeA);

x = A(1,:);
y = A(2,:);
z = A(3,:);
fclose(fid);

fid = fopen(['plotData/BeTSSi_M1_HWBC_BI/FOI_f' num2str(f*100) '_240.dat'], 'wt+','b');
fprintf(fid,'theta\t\t\tTS\n');
for alpha_f_Nr = 1:length(x)    
    fprintf(fid,'%1.15f\t%1.15f\n',x(alpha_f_Nr),y(alpha_f_Nr));
end
fclose(fid); 


fid = fopen(['plotData/BeTSSi_M1_HWBC_BI/FOI_f' num2str(f*100) '_300.dat'], 'wt+','b');
fprintf(fid,'theta\t\t\tTS\n');
for alpha_f_Nr = 1:length(x)    
    fprintf(fid,'%1.15f\t%1.15f\n',x(alpha_f_Nr),z(alpha_f_Nr));
end
fclose(fid); 