function A = importDataSet(filename)

fid = fopen(filename);
A = zeros(10000,2);
iRow = 1;
while (~feof(fid)) 
    temp = cell2mat(textscan(fid,'%f %f\n','CommentStyle','%'));
    if isa(temp,'double') && numel(temp) == 2
        A(iRow,:) = temp;
        iRow = iRow + 1;
    end
end
A(iRow:end,:) = [];
fclose(fid);