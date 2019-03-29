close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

%% Test benchmark models
% alpha_s = 240*pi/180;
alpha_s = 0*pi/180;
% beta_s = 30*pi/180;
beta_s = 0;
beta_f = beta_s;
%         beta_f_arr = beta_s;
%         beta_f = 0*pi/180;
%         beta_f = 0;
beta_f_arr = beta_f;
alpha_f = alpha_s;

d_vec = -[cos(beta_s)*cos(alpha_s);
          cos(beta_s)*sin(alpha_s);
          sin(beta_s)];

%         scatteringCase = 'Ray';
scatteringCase = 'BI';
%         scatteringCase = 'Sweep';
switch scatteringCase
    case 'Ray'
        alpha_f_arr = alpha_f;
        beta_f_arr = beta_f;
        f = 10e3; % frequencies in Hertz
    case 'BI'
        delta_alpha = 1;
        alpha_f_arr = (0:delta_alpha:360)*pi/180;
        alpha_f_arr = 0;
        f = 3e3; % frequencies in Hertz
    case 'Sweep'
        alpha_f_arr = alpha_s;
        delta_f = 1;
        Nfreqs = 3;
        c_f = 1524;
        kstart = 0.01;
        kend = 2;
        f = linspace(kstart*c_f/(2*pi),kend*c_f/(2*pi),Nfreqs);
end
model = 'S5';
% model = 'test';

switch model
    case 'S1'
        setS1Parameters
    case 'S3'
        setS3Parameters
    case 'S5'
        setS5Parameters
    case 'S13'
        setS13Parameters
    case 'S15'
        setS15Parameters
    case 'S35'
        setS35Parameters
    case 'S135'
        setS135Parameters
    case 'test'
        P_inc = 1; % Amplitude of incident wave
        rho_s = [7850, 7850, 7850]; % Density of solid
        rho_f = [1000, 1000, 1.2, 1.2]; % Density of fluids
        c_f = [1500, 1500, 340, 340];  % Speed of sound in fluid domains
        t = [0.008, 0.02, 0.05]; % The thickness of the sphere
        E = [210e9, 210e9, 210e9]; % Youngs modulus of elastic material
        nu = [0.3, 0.3, 0.3]; % Poisson ratio of elastic material
        R_o = [5, 1, 0.5];
    case 'IL'
        setIhlenburgParameters

end
SHBC = true;
ESBC = false;
SSBC = false;
R_i = R_o - t;
omega = 2*pi*f; % Angular frequency

defineBCstring

options = struct('d_vec', d_vec,... 
                 'omega', omega, ...
                 'R_i', R_i, ...
                 'R_o', R_o, ...
                 'P_inc', P_inc, ...
                 'E', E, ...
                 'nu', nu, ...
                 'rho_s', rho_s, ...
                 'rho_f', rho_f, ...
                 'c_f', c_f);
             
specialValues = [];
alpha_f_arr = unique(sort([specialValues', linspace(0,2*pi,3000)]));
v = R_o(1)*[cos(beta_f)*cos(alpha_f_arr); cos(beta_f)*sin(alpha_f_arr); sin(beta_f)*ones(size(alpha_f_arr))]';
data = e3Dss(v, options);

k = omega./c_f;
p_inc = @(v) P_inc*exp(1i*k*dot3(v,d_vec));
% plot(180/pi*alpha_f_arr, real(data(1).p+p_inc(v)), 180/pi*alpha_f_arr, 2*real(p_inc(v)))
plot(180/pi*alpha_f_arr, real(data(1).p+p_inc(v)), 180/pi*alpha_f_arr, real(p_inc_(v,d_vec,k,P_inc,R_o)))
legend('analytic','approx')
% xlim([60,240])
% figure(2)
% plot(180/pi*alpha_f_arr, imag(data(1).p+p_inc(v)), 180/pi*alpha_f_arr, 2*imag(p_inc(v)))
% xlim([60,240])
% figure(3)
% plot(180/pi*alpha_f_arr, abs(data(1).p+p_inc(v)))
% xlim([60,240])

% legend('TS')
% xlabel('k')

function p_inc = p_inc_(v,d_vec,k,P_inc,R)

indices = dot3(v,-d_vec) > 0;
p_inc = zeros(size(v,1),1);
p_inc(indices) = 2*P_inc*exp(1i*k*dot3(v(indices,:),d_vec));
indices = setdiff(1:size(v,1),find(indices));
if ~isempty(indices)
    phi = pi/2-acos(dot3(v(indices,:)/R,d_vec));
    for n = 0:10
        s = R*(phi+2*pi*n);
        p_inc(indices) = p_inc(indices) + exp(1i*k*s)./(s+0.5);
        s = R*(pi-phi+2*pi*n);
        p_inc(indices) = p_inc(indices) + exp(1i*k*s)./(s+0.5);
    end
end


end