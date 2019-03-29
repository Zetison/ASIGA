

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'M4'; % BeTSSi model 5A and BeTSSi model 5B

method = 'BEM';
if strcmp(method, 'BEM')
    formulation = 'CCBIE';
end

f = [1e3, 3e3];             % Frequency

M = 4;

alpha_s = 30*pi/180;
beta_s = 30*pi/180;

plot3Dgeometry = 0;

degreeElev = 0;

loopParameters = {'M','f'};
% plot3Dgeometry = true;
collectIntoTasks