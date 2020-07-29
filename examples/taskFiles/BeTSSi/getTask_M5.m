

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

% model = {'M5A', 'M5B'}; % BeTSSi model 5A and BeTSSi model 5B
model = {'M5A'}; % BeTSSi model 5A and BeTSSi model 5B

method = 'BEM';
formulation = 'CCBIE';

% f = [1e3 3e3];             % Frequency
f = 1e3;             % Frequency

M = 2;

alpha_s = 60*pi/180;
beta_s = 0*pi/180;

degree = 2;
plotResultsInParaview = 0;	% Only if scatteringCase == 'Bi'

solveForPtot = true;
prePlot.plot3Dgeometry = 0;
loopParameters = {'M','model','degree'};
collectIntoTasks
% 
% scatteringCase = {'MS'}; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
% collectIntoTasks
