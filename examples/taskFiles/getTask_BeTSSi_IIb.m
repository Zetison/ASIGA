

scatteringCase = {'BI', 'MS'}; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = {'M5A', 'M5B'}; % BeTSSi model 5A and BeTSSi model 5B

method = {'BEM'};
if strcmp(method, 'BEM')
    formulation = {'CCBIE'};
end

f_arr = [1e3 3e3 10e3 30e3];             % Frequency

M = 1:4;

alpha_s = 60*pi/180;
beta_s = 0*pi/180;

degreeElevArr = 2;

solveForPtot = true;
% collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = {'M4'};

M = 4:7;

alpha_s = 30*pi/180;
beta_s = 30*pi/180;

degreeElevArr = 2;
beta_f_arr = beta_s;

collectIntoTasks