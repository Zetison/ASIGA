function convertSTLtoBDF(filename)
fv = stlread([filename '.stl']);

F = fv.faces;
X = single(fv.vertices);
nX = size(X,1);
nF = size(F,1);
T = table(repmat({'GRID'},nX,1),...
           (1:nX).', ...
           zeros(nX,1), ...
           X(:,1), ...
           X(:,2), ...
           X(:,3));
filename2 = [filename '_temp.bdf'];
writetable(T,filename2,'Delimiter',',','QuoteStrings',false,'FileType','text','WriteVariableNames',false)
T = table(repmat({'CTRIA3'},nF,1),...
           (1:nF).', ...
           zeros(nF,1), ...
           F(:,1), ...
           F(:,2), ...
           F(:,3));
filename = [filename '.bdf'];
writetable(T,filename,'Delimiter',',','QuoteStrings',false,'FileType','text','WriteVariableNames',false)

S = fileread(filename);
S2 = fileread(filename2);
delete(filename2)

S = ['$ Created by Jon Vegard Venaas (NTNU)', newline, S2, S];
FID = fopen(filename, 'w');
if FID == -1
    error('Cannot open file %s', filename); 
end
fwrite(FID, S, 'char');
fclose(FID);


fid = fopen(filename, 'a+');
fprintf(fid, 'ENDDATA');
fclose(fid);