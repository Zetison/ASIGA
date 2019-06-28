clear all

alpha_s = 240*pi/180;
beta_s = 30*pi/180;
beta_f = beta_s;
beta_f_arr = beta_s;
d_vec = -[cos(beta_s)*cos(alpha_s);
          cos(beta_s)*sin(alpha_s);
          sin(beta_s)];
delta_alpha = 0.1;
f = 30e3; % frequencies in Hertz
setS5Parameters
R_i = R_o - t;
omega = 2*pi*f; % Angular frequency
SHBC = 1;
SSBC = 0;
ESBC = 0;

defineBCstring

options = struct('d_vec', d_vec, ...
                 'omega', omega, ...
                 'R_i', R_i, ...
                 'R_o', R_o, ...
                 'P_inc', P_inc, ...
                 'E', E, ...
                 'nu', nu, ...
                 'rho_s', rho_s, ...
                 'rho_f', rho_f, ...
                 'c_f', c_f, ...
                 'calc_farField', 1);

delta_alpha = 0.001;
alpha_f_arr = (0:delta_alpha:360)*pi/180;
pathToResults = '../../../results/e3Dss/';

TS_ = @(alpha) computeTS(alpha,beta_f,options);
TS = TS_(alpha_f_arr);
meanTS = mean(TS);
indices = find((TS(1:end-1)-meanTS).*(TS(2:end)-meanTS) < 0);
specialValues = zeros(numel(indices)-1,2);
% options = optimset('Display','iter');
parfor i = 1:(numel(indices)/2-1)
    tic
    alphaIntr1 = alpha_f_arr(indices(2*i-1:2*i));
    alphaIntr2 = alpha_f_arr(indices(2*i:2*i+1));
    specialValues(i,:) = [fminsearchbnd(@(alpha) -TS_(alpha),mean(alphaIntr1),alphaIntr1(1),alphaIntr1(end)), ...
                   fminsearchbnd(@(alpha) TS_(alpha),mean(alphaIntr2),alphaIntr2(1),alphaIntr2(end))];
    toc
end
specialValues = sort(specialValues(:));
alpha_f_arr = sort([specialValues.', alpha_f_arr]);
plot(alpha_f_arr,TS_(alpha_f_arr))
save([pathToResults 'S5_SHBC_F30_specialValues'],'specialValues')


function TS = computeTS(alpha_f_arr,beta_f,options)

v = [cos(beta_f)*cos(alpha_f_arr); cos(beta_f)*sin(alpha_f_arr); sin(beta_f)*ones(size(alpha_f_arr))]'; 
data = e3Dss(v, options);
TS = 20*log10(abs(data.p));

end