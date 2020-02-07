scatteringCase = 'Sweep';

% model = {'PS'};  % Spherical shell
model = 'SS';  % Spherical shell

method = 'BEM';
formulation = 'CBM';

f = [1e1,5e3];             % Wave number for outer fluid domain
% Eigenfrequencies from Zheng2015itb (note that they are multiplied by two, to
% compensate for the half size of the sphere)
eigenValues = 2*[3.141592653589794 % pi
                   6.283185307179587 % 2*pi
                   9.424777960769379 % 3*pi
                   4.493409457909065
                   7.725251836937707
                   5.763459196894550
                   9.095011330476353
                   6.987932000500519
                   8.182561452571243
                   9.355812111042747]';
eigenValues = eigenValues*1500/(2*pi); % convert to frequency
f = sort([eigenValues, f]);

R_o = 0.5; % Outer radius of shell
M = 2;
alpha_s = pi;
beta_s = 0;    


plotFarField = 0;
useStandardErrorNorm  = 0;
calculateSurfaceError = 1;	% Only for spherical shell and if scatteringCase == 'Bi'
plot3Dgeometry = 0;
calculateFarFieldPattern = 0;

% collectIntoTasks
solveForPtot = true;
