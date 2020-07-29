function printWindRoseData(filename)
% printWindRoseData('../LaTeX/tikzExample/data/421402121491400.dat')
[filepath,name] = fileparts(filename);
if true
%     T = readtable(filename,'Format','%{mm/dd/yyyy}D %{HH:mm}D %f %f','FileType','text', 'HeaderLines',11,'ReadVariableNames',false);
%     temp = T.Var4;
%     T.Var4 = T.Var3;
%     T.Var3 = temp;
%     writetable(T,[filename(1:end-4) '.csv'],'Delimiter', ',')
    
    T = readtable(filename,'Format','%{mm/dd/yyyy}D %{HH:mm}D %f %f','FileType','text', 'HeaderLines',11,'ReadVariableNames',false);
    windSpeeds = T.Var3;
    dirs = T.Var4;
else
    T = readtable(filename,'Format','%d %{yyyy-mm-dd HH:mm}D %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','Delimiter', ',', 'HeaderLines',1,'ReadVariableNames',false);
    dirs = T.Var3;
    windSpeeds = T.Var4;
end
% printResultsToFile2([filepath '/' name '_test.txt'], dirs, windSpeeds)
indices = windSpeeds < 0;
windSpeeds(indices) = [];
dirs(indices) = [];
% windSpeeds(indices) = 0;
% dirs(indices) = 0;
A = zeros(16,7);
windTypes = [0, 0.5, 2, 4, 6, 8, 10, inf];
for i = 1:numel(windSpeeds)
    dirsIdx = ceil(mod(dirs(i)+11.25,360)/22.5);
    windIdx = 1;
    while windSpeeds(i) >= windTypes(windIdx)
        windIdx = windIdx + 1;
    end
    windIdx = max(windIdx - 1,1);
%     windIdx = find(and(windTypes(1:end-1) < windSpeeds(i), windSpeeds(i) <= windTypes(2:end)));
        
    A(dirsIdx,windIdx) = A(dirsIdx,windIdx) + 1;
%     fprintf('%d %d\n',dirsIdx-1,windIdx-1)
%     if i > 20
%         keyboard
%     end
end
A = A/numel(windSpeeds)*100;

meanWind = mean(windSpeeds);
[peakFreq, peakDir] = max(sum(A,2));
dir = {'N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW'};
direction = dir{peakDir};
percentCalm = sum(A(:,1));

fid = fopen([filepath '/' name '_formatted.txt'], 'w+','b');
fprintf(fid, 'South Shore Met Station [SSHR MET] (421402121491400)\n');
fprintf(fid, 'Data from U.S. Geological Survey\n');
fprintf(fid, 'Jan-01-2006\n');
fprintf(fid, 'Dec-31-2010\n');
fprintf(fid, 'https://or.water.usgs.gov/cgi-bin/grapher/grapher.pl\n');
fprintf(fid, '%0.3g \n', meanWind);
fprintf(fid, '%0.4g \n', peakFreq);
fprintf(fid, '%s \n', direction);
fprintf(fid, '%0.3g \n', percentCalm);
fprintf(fid, ['%s ' repmat('%12.2g ',1,size(windTypes(1:end-1),2)) '\n'], '  Dir', windTypes(1:end-1));


for i = 1:size(A,1)
    if i == size(A,1)
        ending = '';
    else
        ending = '\n';
    end
    fprintf(fid, ['%5.1f ' repmat('%12.8g ',1,size(A,2)) ending], 22.5*(i-1), A(i,:));
end
fclose(fid); 


