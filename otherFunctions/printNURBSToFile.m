function printNURBSToFile(nurbs, filename)
% fileID = fopen(filename,'w');
% fprintf(fileID, 'NURBS type: %s\n', nurbs.type);
% fclose(fileID);

dlmwrite(filename,nurbs.number,'precision',16)
dlmwrite(filename,nurbs.degree,'-append','precision',16)

for i = 1:numel(nurbs.knots)
    dlmwrite(filename,nurbs.knots{i},'-append','precision',16)
end


coeffs = nurbs.coeffs;

sz = size(coeffs);
M = zeros(sz(2)*sz(3),4);
for i = 1:sz(2)
    for j = 1:sz(3)
        index = i + (j-1)*sz(2);
        M(index,:) = coeffs(:,i,j);
    end
end
dlmwrite(filename,M,'-append','precision',16)