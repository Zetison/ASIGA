
uiopen('C:\Users\Zetison\Dropbox\work\matlab\results\_studies\studies303\_TAP_TSVSalpha.fig',1)
for study_i = 1:numel(studies)  
    study = studies(study_i);
    options = struct('xname',           'alpha',  ...
                     'yname',           'TS', ...
                     'plotResults', 	1, ... 
                     'printResults',	0, ... 
                     'axisType',        'plot', ... 
                     'lineStyle',       '-', ... 
                     'legendEntries',   {{'method','coreMethod','M','degree'}}, ...
                     'noXLoopPrms',     0); 

    options.xScale = 180/pi;
    figure(1)
    printResultsToTextFiles(study,options)
    
%     study = studies(study_i);
%     options = struct('xname',           'parm',  ...
%                      'yname',           'cond_number', ...
%                      'plotResults', 	1, ... 
%                      'printResults',	0, ... 
%                      'axisType',        'semilogy', ... 
%                      'xLoopName',       'parm', ...
%                      'lineStyle',       '-', ... 
%                      'legendEntries',   {{'method','coreMethod','M','degree'}}, ...
%                      'noXLoopPrms',     1); 
% 
%     options.xScale = 180/pi;
%     figure(4)
%     printResultsToTextFiles(study,options)
end
return
% close all
alpha = (0:0.1:180)*pi/180;
theta_i = alpha;
theta_r = alpha;
setTAPParameters
% R_o = 5.5/1.0936;
% L = 120/1.0936;
f = 299;
c = c_f(1);
omega = 2*pi*f;
k = omega/c;
lambda = 2*pi/k;
TS = TS_BS_submarine(R_o, L, lambda, theta_i, theta_r, E, nu, rho_s, c, false);
TS2 = TS_BS_submarine(R_o, L, lambda, theta_i, theta_r, E, nu, rho_s, c, true);
% sigma_fs = (2*R_o*L/lambda)^2*abs(sin(theta_i).*sin(theta_r)).*bessel_s(0,Alpha,1).^2;

TS_Hodges = importdata('results/TAP/TS.csv');
TS_Hodges = [TS_Hodges; flipud(TS_Hodges(1:end-1,:))];
TS_Hodges(ceil(end/2):end,1) = 180-TS_Hodges(ceil(end/2):end,1);

% % setTAPParameters
% f = 300;
% omega = 2*pi*f;
% k = omega/c;
% lambda = 2*pi/k;
% % theta_i = 0:pi/10000:pi;
% % theta_r = theta_i;
% 
% for j = 1:length(theta_i),
% %     for k = 1:length(th_r),
%         
%         [TS(j)] = TS_BS_submarine(R_o, L, lambda, theta_i(j), theta_r(j));
%         
% %     end
%     
% end

hold on
% plot(theta_i*180/pi, 10*log10(TS),'DisplayName','Karl Thomas') %, alpha*180/pi,TS2)
plot(theta_i*180/pi, TS,'DisplayName','Hodges formler')
plot(theta_i*180/pi, TS2,'DisplayName','Hodges formler uten faktoren 1/2 for endcaps')
% plot(TS_Hodges(:,1), TS_Hodges(:,2),'DisplayName','Rodgers2010 Figure 9.10, page 177') %, alpha*180/pi,TS2)
legend show
% xlabel('Aspect (deg)')
% ylabel('TS')
% 
% figure
% theta_i = (90:0.1:180)*pi/180;
% theta_r = (0:0.1:360)*pi/180;
% [~, TS] = TS_BS_submarine(R_o, L, lambda, theta_i.', theta_r, E, nu, rho_s, c, false);
% imagesc(theta_i*180/pi, theta_r*180/pi, TS.')
% ha = xlabel('Incident angle [deg]');
% set(ha, 'fontweight', 'bold', 'fontsize', 12)
% ha = ylabel('Reflected angle [deg]');
% set(ha, 'fontweight', 'bold', 'fontsize', 12)
% set(gca, 'ydir', 'normal')
% colorbar
% set(gca, 'clim', [-20,30])
% % xlim([90,180])
% % ylim([-180,180])
% box on
