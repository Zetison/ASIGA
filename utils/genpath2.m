
function p = genpath2(d, pattern)
%GENPATH2 calls genpath and removes folders matching a specified pattern
%
% INPUTS:
%   d       ~ char, the name of of the base folder
%   pattern ~ char/cell, the starting folder pattern to exclude
%
% OUTPUTS:
%   p       ~ char, the cleaned path
%
% USAGE:
%   genpath2(folderName) returns an array identical to genpath(folderName)
%   genpath2(folderName, '.git') returns a array without folders starting with .git
%   genpath2(folderName, {'.git', '.svn'}) returns a vector without folders starting with .git or .svn
%
% Santiago I. Sordo-Palacios, 2019
% Call MATLAB's genpath()
p = genpath(d);
% Return if missing or empty input argument
if nargin < 2 || isempty(pattern)
    return
end
% Find folders that match the pattern
splitP = split(p, pathsep);
pattern = strcat(filesep, pattern);
hasPattern = contains(splitP, pattern);
% Index out folders with pattern
cleanP = splitP(~hasPattern);
% Return as list in genpath format
p = char(strjoin(cleanP, pathsep));
end % function-genpath2
