close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

% pathToResults = '../../../results/e3Dss/';
pathToResults = '../results';
set(0,'defaultTextInterpreter','latex');
startMatlabPool
plotP_inc = false;

% applyLoad = 'planeWave';
applyLoad = 'radialPulsation';

f_c = 1500; % (300) source pulse center freq.
N = 2^11;
T = 120/f_c; %60/f_c
B = N/T; % bandwidth
f_L = -B/2;
f_R = B/2;
df = 1/T;
f = linspace(0,f_R-df,N/2);
omega = 2*pi*f;
type = 1;
d_vec = [0,0,1].';
omega_c = 2*pi*f_c;
c_f = 1500;
k_c = omega_c/c_f;

%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotP_inc
    figure(1); clf;
    %             t = linspace(-1/f_c,6/f_c,1000);
    tt = linspace(0,1/f_c,1000);
    tt = [-1e-3,tt,3e-3];
    Pt_inc = Pt_inc_(tt,0,omega_c,k_c,type);
    subplot(311), plot(tt,Pt_inc)
    ft = linspace(0,f_R,2000);
    omegat = 2*pi*ft;
    P_inc = P_inc_(omegat,omega_c,type);
    subplot(312), plot(ft,abs(P_inc))
    subplot(313), plot(ft,atan2(imag(P_inc),real(P_inc)));
    drawnow
    %             printResultsToFile([pathToResults 'Pt_inc'], tt.', Pt_inc.', [], 0, 1)
    %             printResultsToFile([pathToResults 'P_inc'], omegat.', abs(P_inc).', [], 0, 1)
    figure(2)
    plot(omegat,abs(P_inc))
    ylabel('$|P_{\mathrm{inc}}(\omega)|$ [Pa]')
    xlabel('$\omega$ [$\mathrm{s}^{-1}$]')
    savefig([pathToResults 'Figure17b'])

    figure(3)
    plot(tt,Pt_inc)
    ylabel('$\breve{P}_{\mathrm{inc}}(t)$ [Pa]')
    xlabel('$t$ [s]')
    savefig([pathToResults 'Figure17a'])
end
%             return
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot dynamic time domain solution in 1D
SHBC = 1;
SSBC = 0;
ESBC = 0; 
setS15Parameters

R_i = R_o - t;
defineBCstring
npts = 1000;
R_a = 1.5*R_o(1);
%             z = linspace(-1.3*R_a,-R_o(1),npts).';
z1 = linspace(R_o(1),20,npts).';
z2 = linspace(R_i(1),R_o(1),npts).';
z3 = linspace(R_o(2),R_i(1),npts).';
vv1 = [zeros(npts,1), zeros(npts,1), z1];
vv2 = [zeros(npts,1), zeros(npts,1), z2];
vv3 = [zeros(npts,1), zeros(npts,1), z3];
P_inc = @(omega) P_inc_(omega, omega_c,type);

options = struct('d_vec', d_vec, ...
                 'omega', omega(2:end), ...
                 'R_i', R_i, ...
                 'R_o', R_o, ...
                 'P_inc', P_inc, ...
                 'applyLoad',applyLoad,...
                 'calc_sigma_rr',true,...
                 'E', E, ...
                 'nu', nu, ...
                 'rho_s', rho_s, ...
                 'rho_f', rho_f, ...
                 'c_f', c_f);
data = e3Dss({vv1,vv2,vv3}, options);

startIdx = 2000;
startIdx = 2*round(1900/2);
totField1 = zeros(npts,N/2);
totField2 = zeros(npts,N/2);
totField3 = zeros(npts,N/2);
PincField = zeros(npts,N/2);

for n = 0:N-1
    f = f_L + (f_R-f_L)/N*n;
    omega = 2*pi*f;
    k = omega/c_f(1);
    k_vec = d_vec*k;
    if n >= N/2+1
        PincField(:,n-N/2+1) = P_inc_(omega,omega_c,type)*R_o(1).*exp(-1i*k*(norm2(vv1)-R_o(1)))./norm2(vv1);
        totField1(:,n-N/2+1) = PincField(:,n-N/2+1) + data(1).p(:,n-N/2);
        totField2(:,n-N/2+1) = data(1).sigma_rr(:,n-N/2);
        totField3(:,n-N/2+1) = data(2).p(:,n-N/2);
    end
end
dt = T/N;
totFieldTime{1} = 2/T*fft(totField1,N,2);
totFieldTime{2} = 2/T*fft(totField2,N,2);
totFieldTime{3} = 2/T*fft(totField3,N,2);
PincFieldTime = 2/T*fft(PincField,N,2);

for i = 1:3
    temp = totFieldTime{i};
    totFieldTime{i}(:,1:N-startIdx+1) = temp(:,startIdx:end);
    totFieldTime{i}(:,N-startIdx+2:end) = temp(:,1:startIdx-1);
end
temp = PincFieldTime;
PincFieldTime(:,1:N-startIdx+1) = temp(:,startIdx:end);
PincFieldTime(:,N-startIdx+2:end) = temp(:,1:startIdx-1);

figure(2)
m_arr = 0:N-1;
for m = m_arr
    t = dt*m;
    plot(z3,real(totFieldTime{3}(:,m+1)),'blue','DisplayName','p_{tot}');
    hold on
    plot(z2,-real(totFieldTime{2}(:,m+1)),'blue','DisplayName','-\sigma_{rr}');
    plot(z1,real(totFieldTime{1}(:,m+1)),'blue','DisplayName','p_{tot}');
    plot(z1,real(PincFieldTime(:,m+1)),'red','DisplayName','p_{inc}');
    ylim([-2,2])
    plot(R_o(1)*[1,1],ylim,'--','color','black')
    plot(R_i(1)*[1,1],ylim,'--','color','black')
    plot(R_o(2)*[1,1],ylim,'--','color','black')
    hold off
    titleStr = sprintf('Time step %4d, t = %5.3fs. T = %5.3fs.',m,t,T);
    title(titleStr)
    legend show
    drawnow
end
