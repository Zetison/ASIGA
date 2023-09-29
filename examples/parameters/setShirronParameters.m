function varCol = setShirronParameters()

varCol{1}.rho_f = 1000; % Density of fluids
varCol{1}.c_f = 1500;  % Speed of sound in fluid domains
R_o = 1;
L = 10*R_o;
% L = 4*R_o; % Shirron example 1
varCol{1}.R_o = R_o;
varCol{1}.L = L; % Shirron example 2
varCol{1}.mult = round(L/(R_o*pi/2));