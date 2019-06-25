close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

pathToResults = '../../../plotData/e3Dss/';
% pathToResults = '../results';

startMatlabPool

%% Ihlenburg (1998) example
Ieye = [1, 0;
        0, 1;
        0, 0];
ESBC = 0;
for i = 1:3
    if i == 1
        nFreqs = 2000;
        color = [0,70,147]/255;
        legendEntry = 'Sound-hard boundary condition';
    elseif i == 2
        nFreqs = 5000;
        color = [178,0,0]/255;
        legendEntry = 'Sound-soft boundary condition';
    else
        nFreqs = 5000;
        color = [59,124,37]/255;
        legendEntry = 'Neumann-Neumann boundary condition';
    end
    P_inc = 1; % Amplitude of incident wave
    rho_f = [1000, 1000]; % Density of outer fluid
    rho_s = 7669; % Density of solid
    c_f = [1524, 1524];   % Speed of sound in outer (exterior) fluid domain
    t = 0.15; % The thickness of the sphere
    R = 5; % Midsurface radius
    R_o = R+t/2; % Outer radius of shell
    R_i = R-t/2; % Inner radius of shell
    E = 207e9; % Youngs modulus of elastic material
    nu = 0.3; % Poisson ratio of elastic material

    SHBC = Ieye(i,1);
    SSBC = Ieye(i,2);
    defineBCstring

    k = linspace(2/nFreqs,2,nFreqs)'; % wave number
    omega = k*c_f(1);   % Wave number for outer fluid domain

    theta = 180*pi/180;
    d_vec = [1,0,0];
    options = struct('d_vec', d_vec, ...
                     'omega', omega, ...
                     'R_i', R_i, ...
                     'R_o', R_o, ...
                     'P_inc', P_inc, ...
                     'E', E, ...
                     'nu', nu, ...
                     'rho_s', rho_s, ...
                     'rho_f', rho_f, ...
                     'c_f', c_f,...
                     'calc_farField', 1);

    if SHBC
        specialValues = [];
    else
        if false
            v = R_o*[cos(0),0,0];
            f = @(k)-objFunc(k,options,v,c_f(1),0);
            specialValues = findExtremas(f, 2/nFreqs, 2, 100000)';
            v = R_o*[cos(pi),0,0];
            f = @(k)-objFunc(k,options,v,c_f(1),0);
            specialValues = [specialValues; findExtremas(f, 2/nFreqs, 2, 100000)'];
            save([pathToResults 'Ihlenburg_' BC '_extremas'], 'specialValues')
        else
            load([pathToResults 'Ihlenburg_' BC '_extremas'])
        end
        delta = 1e-4;
        specialValues = sort([specialValues; (specialValues-delta); (specialValues+delta)]);
    end
    k = unique(sort([k; specialValues]));
    omega = k*c_f(1);   % Wave number for outer fluid domain
    options.omega = omega;

    v = R_o*[cos(pi),0,0;
             cos(0),0,0];
    data = e3Dss(v, options);

    figure(3)
    F = data(1).p;
    TS = 20*log10(abs(F));
%     plot(k*R_o, TS(1,:))
    plot(k*R_o, TS(1,:),'DisplayName',legendEntry,'color',color)
    set(0,'defaulttextinterpreter','latex')
    hold on
%     title('Ihlenburg (1998) example, $$\theta = 180^\circ$$')
    xlabel('$$k_1 R_{0,1}$$')
    xlim([0, max(k*R_o)])
    ylim([-50, 35])
    ylabel('TS [dB]')  
    legend('off');
    legend('show','location','southeast');
    savefig([pathToResults 'Figure9a'])

    figure(4)
    F = data(1).p;
    TS = 20*log10(abs(F));
    plot(k*R_o, TS(2,:),'DisplayName',legendEntry,'color',color)
    set(0,'defaulttextinterpreter','latex')
    hold on
%     title('Ihlenburg (1998) example - $$\theta = 0^\circ$$')
    xlabel('$$k_1 R_{0,1}$$')
    xlim([0, max(k*R_o)])
    ylim([-50, 35])
    ylabel('TS [dB]')  
    legend('off');
    legend('show','location','southeast');
    savefig([pathToResults 'Figure9b'])
    
    folderName = '../results';
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
%             
%             elevation = '0';
%             frequency = 'S';
%             
%             varCol.alpha_s = pi;
%             varCol.beta_s = 0;
%             scatteringCase = 'Sweep';
%             model = 'IL';
%             varCol.scatteringCase = scatteringCase;
%             
%             varCol.f_arr = omega/(2*pi);
%             
%             saveName = [model '_' BC '_' scatteringCase '_A180_E' elevation '_F' frequency];
%             varCol.saveName = saveName;
%             filename = [folderName '/' saveName];
%             printResultsToFile(filename, k*R_o, TS(1,:).', varCol, 1, 0, 'NTNU_FFI', 'Analytic solution')
% %             
%             saveName = [model '_' BC '_' scatteringCase '_A0_E' elevation '_F' frequency];
%             varCol.saveName = saveName;
%             filename = [folderName '/' saveName];
%             printResultsToFile(filename, k*R_o, TS(2,:).', varCol, 1, 0, 'NTNU_FFI', 'Analytic solution')
% %             
    if 0
        figure(40+i)
        nFreqs = 500;
        k = linspace(2/nFreqs,2,nFreqs)'; % wave number
        k = unique(sort([k; specialValues]));
        omega = k*c_f(1);   % Wave number for outer fluid domain
        options.omega = omega;

        createConvergencePlot('3D',options,v,35, [pathToResults 'IhlenburgError_' num2str(i)])
        savefig([pathToResults 'Figure' num2str(9+i)])
    end
end

% figure
