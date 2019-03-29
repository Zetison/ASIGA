close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

pathToResults = '../../plotData/e3Dss/';
% pathToResults = '../results';

startMatlabPool

%% Chang and Demkowiz (1994) example
P_inc = 1; % Amplitude of incident wave
rho_f = 1000; % Density of fluids
c_f = 1460; % Speed of sound in fluid domains
rho_s = 7800; % Density of solid
E = 2.0e11;
nu = 0.3;
a = 1; % mid surface radius
h = 0.01*a;
R_i = a - h/2; % Inner radius of shell
R_o = a + h/2;  % Outer radius of shell

%%%%%%%%%
k = [15, 20];
omega = k*c_f(1);

d_vec = [0,0,1].';
p_inc = @(v) P_inc*exp(1i*dot3(v,d_vec)*k);
p_tot_Chang15 = importdata('../models/Chang1994voa2/Figure16.csv');
p_tot_Chang20 = importdata('../models/Chang1994voa2/Figure17.csv');
theta = linspace(0,pi,2000).';
phi = 0;
v = R_o*[sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
options = struct('d_vec', d_vec, ...
                 'omega', omega, ...
                 'R_i', R_i, ...
                 'R_o', R_o, ...
                 'P_inc', P_inc, ...
                 'E', E, ...
                 'nu', nu, ...
                 'rho_s', rho_s, ...
                 'rho_f', rho_f, ...
                 'c_f', c_f);

data = e3Dss(v, options);
p_tot = data(1).p + p_inc(v);

figure(16)
plot(theta*180/pi, real(p_tot(:,1)), p_tot_Chang15(:,1), p_tot_Chang15(:,2))
% title(sprintf('$$h/a = %.2f, k = %.1f, N_{eps} = %d$$', h/a, k(1), data(1).N_eps(1)),'interpreter','latex')
xlabel('$\vartheta$','interpreter','latex')
ylabel('Real part of pressure [Pa]')
ylim([-2 2])
legend({'Present work', 'Reference Solution from Chang (1994)'},'Location','northwest')
xlim([0 180])
xtickformat('degrees')
savefig([pathToResults 'Figure7a'])

figure(17)
plot(theta*180/pi, real(p_tot(:,2)), p_tot_Chang20(:,1), p_tot_Chang20(:,2))
% title(sprintf('$$h/a = %.2f, k = %.1f, N_{eps} = %d$$', h/a, k(2), data(1).N_eps(2)),'interpreter','latex')
xlabel('$\vartheta$','interpreter','latex')
xtickformat('degrees')
ylabel('Real part of pressure [Pa]')
ylim([-2 2])
legend({'Present work', 'Reference Solution from Chang (1994)'})
xlim([0 180])
xtickformat('degrees')
savefig([pathToResults 'Figure7b'])


folderName = '../results';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

aspect = '0';
elevation = 'm90';
waveNumber = num2str(15);
BC = 'SSBC';

varCol.alpha_s = 0;
varCol.beta_s = -90;
scatteringCase = 'BI';
model = 'Chang';
varCol.scatteringCase = scatteringCase;


% 
% varCol.f = k(1)*c_f(1)/(2*pi);
% saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_k' waveNumber '_Chang1'];
% varCol.saveName = saveName;
% filename = [folderName '/' saveName];
% printResultsToFile(filename, p_tot_Chang15(:,1), p_tot_Chang15(:,2), varCol, 1, 0, 'Chang', 'Results using WebPlotDigitizer to scan the results from Chang (1994)')
% 
% 
% varCol.f = k(2)*c_f(1)/(2*pi);
% waveNumber = num2str(20);
% saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_k' waveNumber '_Chang2'];
% varCol.saveName = saveName;
% filename = [folderName '/' saveName];
% printResultsToFile(filename, p_tot_Chang20(:,1), p_tot_Chang20(:,2), varCol, 1, 0, 'Chang', 'Results using WebPlotDigitizer to scan the results from Chang (1994)')
% 
% 
% varCol.f = k(1)*c_f(1)/(2*pi);
% waveNumber = num2str(15);
% saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_k' waveNumber '_1'];
% varCol.saveName = saveName;
% filename = [folderName '/' saveName];
% printResultsToFile(filename, theta*180/pi, real(p_tot(:,1)), varCol, 1, 0, 'NTNU_FFI', 'Analytic solution')
% 
% varCol.f = k(2)*c_f(1)/(2*pi);
% waveNumber = num2str(20);
% saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_k' waveNumber '_2'];
% varCol.saveName = saveName;
% filename = [folderName '/' saveName];
% printResultsToFile(filename, theta*180/pi, real(p_tot(:,2)), varCol, 1, 0, 'NTNU_FFI', 'Analytic solution')

saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation];
filename = [folderName '/' saveName];
figure(42)
createConvergencePlot('2D',options,v,60,filename)
savefig([pathToResults 'Figure8'])