% addpath('/usr/local/AdvanpixMCT/4.6.0.13135')
addpath(genpath('IGAfunctions'))
addpath(genpath('NURBS'))
addpath SEM
addpath(genpath('NURBSgeometries'))
addpath NURBSmeshes
addpath(genpath('utils'))
if ~exist('../e3Dss', 'dir')
    error('The e3Dss repository must be downloaded (can be obtained from GitHub: https://github.com/Zetison/e3Dss) and placed at the same level as ASIGA (as ../ASIGA)');
end
addpath(genpath('../e3Dss'))
addpath(genpath('postProcessing'))
addpath(genpath('examples'))
addpath(genpath('integration'))
addpath subRoutines
addpath(genpath('../export_fig'))
if ~exist('../export_fig', 'dir')
    warning('The export_fig repository must be downloaded (can be obtained from GitHub: https://github.com/altmany/export_fig) and placed at the same level as ASIGA (as ../ASIGA)');
end
set(0,'DefaultLegendAutoUpdate','off')

homeDir = expanduser('~');
folderName = [homeDir '/results/ASIGA'];
if ~exist(folderName, 'dir')
    error('The folder in which results should be stored does not exist. Please make such a folder and alter the variable folderName in startup.m accordingly.')
end