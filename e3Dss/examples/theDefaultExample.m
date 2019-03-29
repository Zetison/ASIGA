close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The default example
% When calling the e3Dss function without any extra options, the default
% set of parameters will be chosen. In particular, the rigid scattering of
% a unit sphere impinged by a plane wave traveling along the z-axis is
% simulated at f=1kHz using c_f = 1500m/s as the speed of sound. 
% The result is plotted on the surface as a function of the polar angle.

theta_arr = linspace(0,pi,2000)'; % Set of angles used for plotting
v = [sin(theta_arr), zeros(size(theta_arr)), cos(theta_arr)]; % Evaluate physical location of plotting points

data = e3Dss(v); % Compute solution

% Plot real part of the scattered pressure at the surface of the unit sphere
plot(theta_arr, real(data(1).p))
xlim([0,pi])
legend('Modulus of scattered field')
xlabel('$$\theta$$, polar angle','interpreter','latex')
ylabel('$$\mathrm{real}(p_1)$$, real part of scattered pressure','interpreter','latex')