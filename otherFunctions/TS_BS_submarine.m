function [TS_SH, TS_SS] = TS_BS_submarine(R_o, L, lambda, theta_i, theta_r, E, nu, rho_s, c, alternative)

%
% Beregner den bistatiske målstyrken, TS, til en UVB med lengde L, radius a og der
% lyden har bølgelengde lambda og kommer inn med en vinkel th relativ 
% UVBens lengste akse (90 grader er bredsiden).
%
%   [TS] = TS_cylinder(a, L, th)
%
% Merk at vinklene må være i intervallet [0,2*pi]
%
% Basert på TAP-modellen fra Hodges-boka.
%
k = 2*pi/lambda;
Alpha = k*L*(cos(theta_i) + cos(theta_r))/2; 

sigma_cbs = R_o*L^2/lambda*(sin(theta_i).*sin(theta_r))./(sin(theta_i) + sin(theta_r)).*besselj(0,Alpha).^2;
% sigma_cbs = R_o*L^2/lambda*sin(theta_i).*sin(theta_r)./(sin(theta_i) + sin(theta_r)).*bessel_s(0,Alpha,1).^2;

sigma_fs = (2*R_o*L/lambda)^2*abs(sin(theta_i).*sin(theta_r)).*besselj(0,Alpha).^2;
% sigma_fs = (2*R_o*L/lambda)^2*abs(sin(theta_i).*sin(theta_r)).*bessel_s(0,Alpha,1).^2;
if alternative
    sigma_ecbs = (1/1.1^4)*(k*R_o)^4*R_o^2/4 * (k*R_o<=1.1) + R_o^2/4 * (k*R_o>1.1);
else
    sigma_ecbs = (1/1.1^4)/2*(k*R_o)^4*R_o^2/4 * (k*R_o<=1.1) + 1/2*R_o^2/4 * (k*R_o>1.1);
end



indices = abs(acos(cos(theta_i).*cos(theta_r)+sin(theta_i).*sin(theta_r))) <= pi/2;
% return
sigma_s = sigma_cbs .* indices + sigma_fs .* ~indices; % Tolker det dithen at en enten har BS eller FS, ikke begge deler samtidig.

sigma_s(abs(sin(theta_i).*sin(theta_r)) < 1e3*eps) = 0;

nu2 = 2; % skaleringsfaktor
K = E/(3*(1-2*nu));
G = E/(2*(1+nu));
c_nu = sqrt((3*K+4*G)/(3*rho_s));

% c_nu = 5667.814911773838; % Lydhastigheten i skroget (Hodges er litt vinglete på om det er trykkhastigheten eller skjærhastigheten som skal brukes her. Jeg går for skjær (3000m/s), trykkhastigheten er 6000m/s.
B = 0.2; % effektivitetsfaktor på typisk 0.2
dth_i = abs(theta_i-pi) - pi/2;
dth_r = abs(theta_r-pi) - pi/2;
th_lim = acos(c/c_nu);

P_i = B*L/pi*cos(pi/(2*nu2)*dth_i./(pi/2-th_lim));
P_i(abs(dth_i) >= nu2*abs(pi/2 - th_lim)) = 0;
P_r = B*L/pi*cos(pi/(2*nu2)*dth_r./(pi/2-th_lim));
P_r(abs(dth_r) >= nu2*abs(pi/2 - th_lim)) = 0;

sigma_ew = P_i.*P_r;

TS_SH = 10*log10(sigma_s + sigma_ecbs);
TS_SS = 10*log10(abs(sigma_s + sigma_ecbs + sigma_ew));

% if (1),
%     figure; plot(th_i, 10*log10(abs([sigma_cbs;sigma_ecbs;sigma_fs;sigma_ew])))
% end
    
return