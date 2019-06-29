

scatteringCase = {'BI'}; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = {'M5A', 'M5B'}; % BeTSSi model 5A and BeTSSi model 5B

method = {'BEM'};
if strcmp(method, 'BEM')
    formulation = {'CCBIE'};
end

f_arr = [1e3 3e3];             % Frequency

M = 1;

alpha_s = 60*pi/180;
beta_s = 0*pi/180;

degreeElevArr = 2;
plotResultsInParaview = 1;	% Only if scatteringCase == 'Bi'

solveForPtot = true;
% plot3Dgeometry = true;
collectIntoTasks
% 
% scatteringCase = {'MS'}; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
% collectIntoTasks
