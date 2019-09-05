function generate_g2_file(nurbs,outputFileName)
patch = 1;
parmDim = numel(nurbs{patch}.knots);
controlPts = nurbs{patch}.coeffs;  
d = size(controlPts,1)-1;
switch parmDim
    case 1
        classType = 100;
    case 2
        classType = 200;
    case 3
        classType = 700;
end
    	

fid = fopen(outputFileName,'wt+','b');
fprintf(fid,'%d 1 0 0\n',classType);
fprintf(fid,'%d 1\n', d);    % latter number: 1=rational, 0=non-rational



for i = 1:parmDim
    n = nurbs{patch}.number(i);
    p = nurbs{patch}.degree(i);
    Xi = nurbs{patch}.knots{i};
    fprintf(fid,'%d %d\n',n,p+1);
    for j = 1:length(Xi)
        fprintf(fid, '%g ',Xi(j));
    end
    fprintf(fid, '\n');
end
 
for i = 1:d
    controlPts(i,:,:,:) = controlPts(i,:,:,:).*controlPts(end,:,:,:);
end
controlPts = reshape(controlPts,d+1,[]);

for i = 1:size(controlPts,2)
    fprintf(fid, [repmat('%20.15g',1,size(controlPts,1)) '\n'], controlPts(:,i));
end
fclose(fid);

