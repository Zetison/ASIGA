function T = strctCp(S, fieldList)
if nargin == 1
   fieldList = fieldnames(S);
end 
for iField = 1:numel(fieldList)
   Field    = fieldList{iField};
   T.(Field) = S.(Field);
end