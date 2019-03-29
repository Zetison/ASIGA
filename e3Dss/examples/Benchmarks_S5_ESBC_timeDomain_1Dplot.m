close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

pathToResults = '../../plotData/e3Dss/';
% pathToResults = '../results';

startMatlabPool

if 0
    f_c = 300; % (300) source pulse center freq.
    N = 2^9;
    M = 2*N; % corresponds to fs = sampling rate
    T = 30/f_c;
    B = N/T; % bandwidth
    f_L = -B/2;
    f_R = B/2;
    df = 1/T;
    type = 1;
    useNegativeFreqsInScatPres = 0;
    if useNegativeFreqsInScatPres
        f = linspace(f_L,f_R-df,N);
    else
        f = linspace(0,f_R-df,N/2);
    end
    omega = 2*pi*f;

    d_vec = [0,0,1].';
    R_o = 2;
    omega_c = 2*pi*f_c;
    P_inc = @(omega) P_inc_(omega, omega_c,type);
    c_f = 1500;
    k_c = omega_c/c_f;
    npts = 400;
    z = linspace(-40,-R_o,npts).';

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plots
    figure(1); clf;
    t = linspace(0,30/f_c,1000);
    subplot(311), plot(t,Pt_inc_(t,0,omega_c,k_c,type))
    Sw1 = P_inc_(omega,omega_c,type);
    subplot(312), plot(f,abs(Sw1))
    subplot(313), plot(f,atan2(imag(Sw1),real(Sw1)));
    drawnow
%         return
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot dynamic time domain solution in 1D
    vv = [zeros(npts,1), zeros(npts,1), z];
    options = struct('d_vec', d_vec, ...
                     'omega', omega(2:end), ...
                     'R_o', R_o, ...
                     'P_inc', P_inc, ...
                     'c_f', c_f);
    data = e3Dss(vv, options);
else
    f_c = 1500; % (300) source pulse center freq.
    N = 2^10;
    M = N; % corresponds to fs = sampling rate
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
    % Plots
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
            return
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot dynamic time domain solution in 1D
%             setS15Parameters
    setS5Parameters
    npts = 1000;
    R_a = 1.5*R_o(1);
%             z = linspace(-1.3*R_a,-R_o(1),npts).';
    z = linspace(R_o(1),20,npts).';
    R_i = R_o - t;
    ESBC = 1;
    vv = [zeros(npts,1), zeros(npts,1), z];
    rho_f = rho_f(1:end-1);
    c_f = c_f(1:end-1);
    R_i = R_i(1:end-1);
    P_inc = @(omega) P_inc_(omega, omega_c,type);

    options = struct('d_vec', d_vec, ...
                     'omega', omega(2:end), ...
                     'R_i', R_i, ...
                     'R_o', R_o, ...
                     'P_inc', P_inc, ...
                     'E', E, ...
                     'nu', nu, ...
                     'rho_s', rho_s, ...
                     'rho_f', rho_f, ...
                     'c_f', c_f);
    data = e3Dss(vv, options);
end
startIdx = 2000;
totField = zeros(npts,N/2);
PincField = zeros(npts,N/2);

for n = 0:N-1
    f = f_L + (f_R-f_L)/N*n;
    omega = 2*pi*f;
    k = omega/c_f(1);
    k_vec = d_vec*k;
    if n >= N/2+1
        PincField(:,n-N/2+1) = P_inc_(omega,omega_c,type).*exp(1i*dot3(vv, k_vec));
        totField(:,n-N/2+1) = PincField(:,n-N/2+1) + data(1).p(:,n-N/2);
    end
end
dt = T/M;
totFieldTime = 2/T*fft(totField,M,2);
PincFieldTime = 2/T*fft(PincField,M,2);

temp = totFieldTime;
totFieldTime(:,1:M-startIdx+1) = temp(:,startIdx:end);
totFieldTime(:,M-startIdx+2:end) = temp(:,1:startIdx-1);
temp = PincFieldTime;
PincFieldTime(:,1:M-startIdx+1) = temp(:,startIdx:end);
PincFieldTime(:,M-startIdx+2:end) = temp(:,1:startIdx-1);

figure(2)
m_arr = 0:M-1;
for m = m_arr
    t = dt*m;
    plot(z,real(totFieldTime(:,m+1)),z,real(PincFieldTime(:,m+1)));
    titleStr = sprintf('Time step %4d, t = %5.3fs. T = %5.3fs.',m,t,T);
    title(titleStr)
    ylim([-1.3,1.3])
    legend('p_{tot}', 'p_{inc}')
    drawnow
    pause(0.1)
    if m == 120
        pause
    end
end
figure(3)
plot(dt*m_arr,real(totFieldTime(1,:)),dt*m_arr,real(PincFieldTime(1,:)));
xlim([0, 120*dt])
