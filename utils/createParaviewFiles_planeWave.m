close all
clear all

addpath IGAfunctions
addpath NURBS
addpath NURBSgeometries
addpath postProcessing
addpath otherFunctions
addpath e3Dss
addpath e3Dss/utils
addpath e3Dss/models
addpath integration
addpath subRoutines
addpath export_fig
startMatlabPool


f = 1000;
f_c = 1000;
L = 5; %5, 0.5
N = 2^9;
omega = 2*pi*f;
P_inc = 1;
air = 1;
M = 160;
if air
    h = 0.02; % thickness of shell
    rhof2 = 1.2;
    cf2 = 340;
else
    h = 0.08; % thickness of shell
    rhof2 = 1000;
    cf2 = 1500;
end
% h = L/2; % thickness of shell
E = 210e9;
nu = 0.3;
rhos = 7850;
rhof1 = 1000;
cf1 = 1500;

theta1 = 10*pi/180;
options = struct('h',   	h,  ...      % Thickness of solid layer
                 'omega',   omega, ...   % Angular frequency
                 'theta1', 	theta1,   ...     	% Incident angle
                 'rhof1', 	rhof1,    ...    % Mass density of upper fluid domain
                 'rhof2', 	rhof2,    ...    	% Mass density of lower fluid domain
                 'P_inc',   P_inc,    ...     	% Amplitude of incident wave
                 'E',       E,   ...    % Youngs modulus for solid layers
                 'nu',      nu,   ...    	% Poisson ratio for solid layers
                 'cf2',     cf2,   ...    	% Speed of sound in lower fluid domain
                 'rhos',    rhos,   ...     % Mass densities of solid layer
                 'cf1',     cf1);          % Speed of sound in upper fluid domain

[visElements{1}, nodes{1}] = meshRectangle([-L, 0],[L, L/2-h/2],M);
[visElements{2}, nodes{2}] = meshRectangle([-L, -h],[L, 0],M);
[visElements{3}, nodes{3}] = meshRectangle([-L, -h/2-L/2],[L, -h],M);
d_vec = [sin(theta1), -cos(theta1)]';
folderName = 'plotData/planeScattering/paraviewResults/';

c_f = [cf1,cf2];
rho_f = [rhof1,rhof2];
N_fine = 2*N;
% T = (L/2-h/2)/c_f(1)+h/min(cs1,cs2)+(L/2-h/2)/c_f(2)+L/c_f(1); %N/M;
T = 0.008; %N/M;
B = N/T; % bandwidth

f_R = B/2;
df = 1/T;
f = linspace(0,f_R-df,N/2);
omega = 2*pi*f;
options.omega = omega(2:end);
startIdx = 640; % 900
% startIdx = 1; % 900

K = E/(3*(1-2*nu));
G = E/(2*(1+nu));
cs1 = sqrt((3*K+4*G)/3/rhos);
cs2 = sqrt(G/rhos);
npts = 0;
for i = 1:length(nodes)
    npts = npts + size(nodes{i},1);
end
type = 1;
totnpts = npts*N*4
omega_c = 2*pi*f_c;
options.P_inc = @(omega) P_inc_(omega,omega_c,type);

data = planeScattering(nodes,options);

% omega = 1;
% options.omega = omega;
% X = [pi,0];
% Y = [pi,-h];
% k = omega/c_f(1);
% k_vec = d_vec*k;
% data = planeScattering({X,[X;Y],Y},options);
% abs((data(1).dPdz(1) + 1i*P_inc_(omega,omega_c,type)*k_vec(2).*exp(1i*dot3(X, k_vec)))/rho_f(1)/omega^2 - data(1).u_z(1))/abs(data(1).u_z(1))
% abs(data(2).dPdz(1)/rho_f(2)/omega^2 - data(1).u_z(2))/abs(data(1).u_z(2))


[pathstr,filename] = fileparts(folderName);

m = 1;
for j = 1:3
    VTKdata.displacement = zeros(size(nodes{j},1),3,length(omega));
    VTKdata.totField = zeros(size(nodes{j},1),1,length(omega));
    if mod(j,2) == 0
        VTKoptions = struct('name',[pathstr '/solid' num2str(m) filename], 'celltype', 'VTK_TRIANGLE', 'plotTimeOscillation', 1, ...
                            'plotTotField',1,'plotDisplacementVectors',1,'plotVonMisesStress',1); 
        VTKdata.stress = zeros(size(nodes{j},1),6,length(omega));
        for i = 2:length(omega)  
%             VTKdata.totField(:,:,i) = data(m).phi(:,i-1);
            VTKdata.displacement(:,1,i) = data(m).u_x(:,i-1);
            VTKdata.displacement(:,2,i) = data(m).u_z(:,i-1);
            VTKdata.stress(:,1,i) = data(m).sigma_xx(:,i-1);
            VTKdata.stress(:,2,i) = data(m).sigma_yy(:,i-1);
            VTKdata.stress(:,3,i) = data(m).sigma_zz(:,i-1);
            VTKdata.stress(:,5,i) = data(m).sigma_xz(:,i-1);
            VTKdata.totField(:,:,i) = data(m).sigma_zz(:,i-1);
%             VTKdata.totField(:,:,i) = data(m).phi(:,i-1) + data(m).psi(:,i-1);
        end
        VTKdata.stress = 2/T*real(fft(VTKdata.stress,N_fine,3));
        tempStress = VTKdata.stress;
        VTKdata.stress(:,:,1:N_fine-startIdx+1) = tempStress(:,:,startIdx:end);
        VTKdata.stress(:,:,N_fine-startIdx+2:end) = tempStress(:,:,1:startIdx-1);
    else
        VTKoptions = struct('name',[pathstr '/fluid' num2str(m) filename], 'celltype', 'VTK_TRIANGLE', 'plotTimeOscillation', 1, ...
                                     'plotSphericalRadialDisplacement',0, 'plotTotField', 1,'plotDisplacementVectors',1);
        if m == 1
            for i = 2:length(omega)  
                k = omega(i)/c_f(m);
                k_vec = d_vec*k;
                p_inc = @(v) P_inc_(omega(i),omega_c,type).*exp(1i*dot3(v, k_vec));
                dp_incdx = @(v) 1i*P_inc_(omega(i),omega_c,type)*k_vec(1).*exp(1i*dot3(v, k_vec));
                dp_incdz = @(v) 1i*P_inc_(omega(i),omega_c,type)*k_vec(2).*exp(1i*dot3(v, k_vec));
                VTKdata.displacement(:,1,i) = (data(m).dPdx(:,i-1) + dp_incdx(nodes{m}))/rho_f(m)/omega(i)^2;
                VTKdata.displacement(:,2,i) = (data(m).dPdz(:,i-1) + dp_incdz(nodes{m}))/rho_f(m)/omega(i)^2;
                VTKdata.totField(:,:,i) = data(m).P(:,i-1) + p_inc(nodes{m});
            end
        else
            for i = 2:length(omega)  
                VTKdata.totField(:,:,i) = data(m).P(:,i-1);
                VTKdata.displacement(:,1,i) = data(m).dPdx(:,i-1)/rho_f(m)/omega(i)^2;
                VTKdata.displacement(:,2,i) = data(m).dPdz(:,i-1)/rho_f(m)/omega(i)^2;
            end
        end
    end
    VTKdata.totField = 2/T*real(fft(VTKdata.totField,N_fine,3));
    VTKdata.displacement = 2/T*real(fft(VTKdata.displacement,N_fine,3));

    if mod(j,2) == 0 
        m = m + 1;
    end
    temp = VTKdata.totField;
    tempDispl = VTKdata.displacement;
    VTKdata.totField(:,:,1:N_fine-startIdx+1) = temp(:,:,startIdx:end);
    VTKdata.totField(:,:,N_fine-startIdx+2:end) = temp(:,:,1:startIdx-1);
    VTKdata.displacement(:,:,1:N_fine-startIdx+1) = tempDispl(:,:,startIdx:end);
    VTKdata.displacement(:,:,N_fine-startIdx+2:end) = tempDispl(:,:,1:startIdx-1);
    VTKdata.nodes = [nodes{j}, ones(size(nodes{j},1),1)];
    VTKdata.visElements = visElements{j};
    VTKdata.omega = omega;
    VTKoptions.T = T;

    VTKoptions.N = N_fine;
    % disp(['Time spent on assembling data ' num2str(toc) ' seconds.'])
    tic
    makeVTKfile(VTKdata, VTKoptions);
    disp(['Time spent on making VTK file ' num2str(toc) ' seconds.'])
end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rho_f = options.rho_f;
% c_f = options.c_f;
% SSBC = options.SSBC;
% noDomains = 2*M+1-ESBC-2*SHBC-SSBC;
% nodes1 = cell(noDomains,1);
% visElements = cell(noDomains,1);
% % keyboard
% % nodes(2:2:end) = zeros(1,0);
% npts = 0;
% for i = 1:length(nodes1)
%     npts = npts + size(nodes1{i},1);
% end
% npts
% % keyboard
% % nodes = nodes(1);
% disp(['Time spent building mesh ' num2str(toc) ' seconds.'])
% % keyboard
% % nodes = nodes(1:end-1);
% % nodes = nodes(1);
% if plotTimeOscillation
%     N = 30;
%     N_fine = 30;
%     type = 1;
%     startIdx = 1;
%     options.P_inc = 1;
% elseif plotInTimeDomain
%     f_c = options.f_c;
%     N = options.N; % 200
%     N_fine = 2*N;
%     T = 2*60/f_c; %N/M;
%     B = N/T; % bandwidth
%     
%     f_R = B/2;
%     df = 1/T;
%     f = linspace(0,f_R-df,N/2);
%     omega = 2*pi*f;
%     options.omega = omega(2:end);
%     if options.usePlaneWave
%         startIdx = 1900; % 900
%     elseif options.usePointChargeWave
%         startIdx = 2000;
%     end
%     type = 1;
%     totnpts = npts*N*4
%     omega_c = 2*pi*f_c;
%     options.P_inc = @(omega) P_inc_(omega,omega_c,type);
% %     keyboard
% else
%     N = 1;
%     N_fine = 1;
%     type = 1;
%     startIdx = 1;
%     options.P_inc = 1;
% end
% tic
% data = e3Dss(nodes1, options);
% disp(['Time spent on computing solution: ' num2str(toc) ' seconds.'])
% m = 1;
% [pathstr,filename] = fileparts(vtfFileName);
% 
% for j = 1:length(nodes1)
%     if ~(SSBC && j == length(nodes1))
%         tic
%         clear VTKdata
%         if mod(j,2) == 0 && computeForSolidDomain
%             toggleJacobianMatrix = false(1,9);
%             if ESBC && m == M
%                 VTKoptions = struct('name',[pathstr '/solid' num2str(m) filename], 'celltype', 'VTK_TRIANGLE', 'plotTimeOscillation', plotTimeOscillation, ...
%                             'plotSphericalRadialDisplacement',1, 'plotdu_xdx',1, 'plotdu_xdy',1, 'plotdu_xdz',1, 'plotdu_ydx',1, 'plotdu_ydy',1, 'plotdu_ydz',1, 'plotdu_zdx',1, 'plotdu_zdy',1, 'plotdu_zdz',1, ...
%                             'plotDisplacementVectors',0,'plotSphericalStress_rr',1,'plotStressXX',1,'plotStressZZ',1,'plotVonMisesStress',1); 
%                 toggleJacobianMatrix = [VTKoptions.plotdu_xdx VTKoptions.plotdu_xdy VTKoptions.plotdu_xdz ...
%                                         VTKoptions.plotdu_ydx VTKoptions.plotdu_ydy VTKoptions.plotdu_ydz ...
%                                         VTKoptions.plotdu_zdx VTKoptions.plotdu_zdy VTKoptions.plotdu_zdz];
%             else
%                 VTKoptions = struct('name',[pathstr '/solid' num2str(m) filename], 'celltype', 'VTK_WEDGE', 'plotTimeOscillation', plotTimeOscillation, ...
%                             'plotSphericalRadialDisplacement',0, 'plotDisplacementVectors',0,'plotSphericalStress_rr',1,'plotVonMisesStress',0); 
%             end
%                     
%             if VTKoptions.plotDisplacementVectors || VTKoptions.plotSphericalRadialDisplacement
%                 displacement = zeros(size(nodes1{j},1),3,length(omega));
%                 for i = 2:length(omega)
%                     displacement(:,:,i) = [data(m).u_x(:,i-1) data(m).u_y(:,i-1) data(m).u_z(:,i-1)];
%                 end
%                 if ~plotInTimeDomain
%                     VTKdata.displacement = [data(m).u_x data(m).u_y data(m).u_z];
%                     
%                     temp = VTKdata.displacement;
%                     VTKdata.displacement = zeros([size(VTKdata.displacement), N]);
%                     for i = 1:N
%                         t = (i-1)/N*2*pi/omega;
%                         VTKdata.displacement(:,:,i) = real(temp*exp(-1i*omega*t));
%                     end
%                 end
%             end
%             if any(toggleJacobianMatrix)
%                 VTKdata.jacobian = zeros(size(nodes1{j},1),6,length(omega));
%                 
%                 if ~plotInTimeDomain
%                     VTKdata.jacobian = [data(m).du_xdx data(m).du_xdy data(m).du_xdz data(m).du_ydx data(m).du_ydy data(m).du_ydz data(m).du_zdx data(m).du_zdy data(m).du_zdz];
%                     
%                     temp = VTKdata.jacobian;
%                     VTKdata.jacobian = zeros([size(VTKdata.jacobian), N]);
%                     for i = 1:N
%                         t = (i-1)/N*2*pi/omega;
%                         VTKdata.jacobian(:,:,i) = real(temp*exp(-1i*omega*t));
%                     end
%                 end
%             end
%             if VTKoptions.plotVonMisesStress || VTKoptions.plotSphericalStress_rr 
%                 VTKdata.stress = zeros(size(nodes1{j},1),6,length(omega));
%                 
%                 if plotInTimeDomain
%                     for i = 2:length(omega)
%                         VTKdata.stress(:,:,i) = [data(m).sigma_xx(:,i-1) data(m).sigma_yy(:,i-1) data(m).sigma_zz(:,i-1) data(m).sigma_yz(:,i-1) data(m).sigma_xz(:,i-1) data(m).sigma_xy(:,i-1)];
%                     end
%                     VTKdata.stress = 2/T*real(fft(VTKdata.stress,N_fine,3));
%                     temp = VTKdata.stress;
%                     VTKdata.stress(:,:,1:N_fine-startIdx+1) = temp(:,:,startIdx:end);
%                     VTKdata.stress(:,:,N_fine-startIdx+2:end) = temp(:,:,1:startIdx-1);
%                 else
%                     VTKdata.stress = [data(m).sigma_xx data(m).sigma_yy data(m).sigma_zz data(m).sigma_yz data(m).sigma_xz data(m).sigma_xy];
%                     
%                     temp = VTKdata.stress;
%                     VTKdata.stress = zeros([size(VTKdata.stress), N]);
%                     for i = 1:N
%                         t = (i-1)/N*2*pi/omega;
%                         VTKdata.stress(:,:,i) = real(temp*exp(-1i*omega*t));
%                     end
%                 end
%             end
%         elseif mod(j,2) ~= 0
%             VTKoptions = struct('name',[pathstr '/fluid' num2str(m) filename], 'celltype', 'VTK_TRIANGLE', 'plotTimeOscillation', plotTimeOscillation, 'plotDisplacementVectors', 0, ...
%                              'plotSphericalRadialDisplacement',0, 'plotTotField', 1, 'plotScalarField', 1, 'plotSPL', 1);
%             if ~(plotInTimeDomain || plotTimeOscillation)
%                 VTKoptions.plotTotFieldAbs = 1;
%             end
%             if VTKoptions.plotDisplacementVectors 
%                 displacement = zeros(size(nodes1{j},1),3,length(omega));
%             end
%             VTKdata.totField = zeros(size(nodes1{j},1),1,length(omega));
%             if m == 1
%                 if plotInTimeDomain
%                     VTKdata.P_inc = zeros(size(nodes1{j},1),1,length(omega));
%                     for i = 2:length(omega)                     
%                         if options.usePlaneWave
%                             k = omega(i)/c_f(1);
%                             k_vec = options.d_vec*k;
%                             p_inc = @(v) P_inc_(omega(i),omega_c,type).*exp(1i*dot3(v, k_vec));
%                         elseif options.usePointChargeWave
%                             r_s = options.r_s;
%                             k = omega(i)/c_f(1);
%                             r = @(y) norm2(repmat(-r_s*d_vec.',size(y,1),1)-y);
%                             Phi_k = @(y) exp(1i*k*r(y))./(4*pi*r(y));     
%                             p_inc = @(y) 4*pi*r_s*P_inc_(omega(i),omega_c,type)*Phi_k(y);
%                         end
%                         VTKdata.P_inc(:,:,i) = p_inc(nodes1{1});
%                         VTKdata.totField(:,:,i) = data(m).p(:,i-1) + VTKdata.P_inc(:,:,i);
%                     end
%                 else
%                     k = omega/c_f(1);
%                     k_vec = options.d_vec*k;
%                     p_inc = @(v) exp(1i*dot3(v, k_vec));   
%                     VTKdata.P_inc = p_inc(nodes1{1});
%                     VTKdata.totField = data(m).p(:,1) + VTKdata.P_inc;
%                     VTKdata.scalarField = data(m).p(:,1);
%                 end
%                 VTKdata = rmfield(VTKdata,'P_inc');
%             else
%                 VTKoptions.plotScalarField = false;
%                 VTKoptions.plotP_inc = false;
%                 if plotInTimeDomain
%                     for i = 2:length(omega)
%     %                     gScalarField = [data(m).dpdx(:,i) data(m).dpdy(:,i) data(m).dpdz(:,i)];
%     %                     displacement(:,:,i) = gScalarField/(rho_f(m)*omega(i)^2);
%                         VTKdata.totField(:,:,i) = data(m).p(:,i-1);
%                     end
%                 else
%                     VTKdata.totField = data(m).p;
%                 end
%             end
%             if plotInTimeDomain
%                 VTKdata.totField = 2/T*real(fft(VTKdata.totField,N_fine,3));
%                 temp = VTKdata.totField;
%                 VTKdata.totField(:,:,1:N_fine-startIdx+1) = temp(:,:,startIdx:end);
%                 VTKdata.totField(:,:,N_fine-startIdx+2:end) = temp(:,:,1:startIdx-1);
%             elseif plotTimeOscillation
%                 temp = VTKdata.totField;
%                 VTKdata.totField = zeros([size(VTKdata.totField), N]);
%                 for i = 1:N
%                     t = (i-1)/N*2*pi/omega;
%                     VTKdata.totField(:,:,i) = real(temp*exp(-1i*omega*t));
%                 end
%             else
%                 VTKdata.totFieldAbs = abs(VTKdata.totField);
%                 VTKdata.totField = real(VTKdata.totField);
%             end
%         end
%         if mod(j,2) == 0 
%             m = m + 1;
%         end
% %         VTKdata.displacement = displacement;
%         VTKdata.nodes = nodes1{j};
%         VTKdata.visElements = visElements{j};
%         VTKdata.omega = omega;
%         if plotInTimeDomain
%             VTKoptions.T = T;
%         end
%         VTKoptions.N = N_fine;
%         disp(['Time spent on assembling data ' num2str(toc) ' seconds.'])
%         tic
%         if mod(j,2) || computeForSolidDomain
%             makeVTKfile(VTKdata, VTKoptions);
%         end
%         disp(['Time spent on making VTK file ' num2str(toc) ' seconds.'])
%     end
% end
%        