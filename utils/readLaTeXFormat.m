function X = readLaTeXFormat(fileName)



fileID = fopen(fileName);
headerLines = {};
i = 1;
X = zeros(10000,2);
while ~feof(fileID)
    line = fgetl(fileID);
    if ~isempty(line) && ~strcmp(line(1),'%')
        if isempty(headerLines)
            headerLines = strsplit(line);
        else
            X(i,:) = str2num(line);
            i = i + 1;
        end
    end
end
X = X(1:i-1,:);