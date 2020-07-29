close all
clear all

addpath export_fig

plotFullSphere = 1;

M = 1;
f_arr = 100;
alpha_s_arr = 240*pi/180;
k_wn_Nr = 1;
load(['plotData/M3_SHBC_BI/FULLf_' num2str(f_arr(k_wn_Nr)) '_mesh' num2str(M) '_' num2str(round(alpha_s_arr*180/pi))])

TS = reshape(TS,size(TS,2),size(TS,3));
if ~plotFullSphere
    indices = find((-20.01 <= beta_f_arr*180/pi).*(beta_f_arr*180/pi <= 20.01));
    beta_f_arr = beta_f_arr(indices);
    TS = TS(:,indices);
end

X = zeros(length(alpha_f_arr),length(beta_f_arr));
Y = zeros(length(alpha_f_arr),length(beta_f_arr));
Z = zeros(length(alpha_f_arr),length(beta_f_arr));

eta2 = 0.8;
eta1 = 0.265;
t = 0.008;   % Thickness
R_o1 = 5; % Outer radius
R_o2 = 3; % Outer radius
R_i1 = R_o1-t; % Inner radius
R_i2 = R_o2-t; % Inner radius
L = 49-R_o1-R_o2;
x_0 = [L/2+(R_o1-R_o2)/2; 0; 0];
nurbs = getModel3Data(R_o1, R_o2, R_i1, R_i2, L, 'Xaxis', eta1, eta2, x_0);

R_far = 1.6*max(max(max(max(abs(nurbs.coeffs(1:3,:,:,:))))));
for alpha_f_Nr = 1:length(alpha_f_arr)
    alpha_f = alpha_f_arr(alpha_f_Nr);
    for beta_f_Nr = 1:length(beta_f_arr)
        beta_f = beta_f_arr(beta_f_Nr);
        
        X(alpha_f_Nr, beta_f_Nr)  = R_far*cos(beta_f)*cos(alpha_f);
        Y(alpha_f_Nr, beta_f_Nr)  = R_far*cos(beta_f)*sin(alpha_f);
        Z(alpha_f_Nr, beta_f_Nr)  = R_far*sin(beta_f);
    end
end   

surf(X,Y,Z,TS,'EdgeColor','none','LineStyle','none')

hold on
if max(beta_f_arr*180/pi) < 30
    plotNURBS(nurbs,[100 100 0], 1, getColor(1), 1);
end
view(-40+90, 30)
colorbar
colormap jet
axis equal
grid off
axis off
hold on  
camlight('left')

if max(beta_f_arr*180/pi) < 30
    alpha_f_arr_meshlines = (0:20:340)*pi/180;
    for alpha_f_Nr = 1:length(alpha_f_arr_meshlines)
        alpha_f = alpha_f_arr_meshlines(alpha_f_Nr);
        x = zeros(size(beta_f_arr));
        y = zeros(size(beta_f_arr));
        z = zeros(size(beta_f_arr));
        for beta_f_Nr = 1:length(beta_f_arr)
            beta_f = beta_f_arr(beta_f_Nr);
            x(beta_f_Nr) = R_far*cos(beta_f)*cos(alpha_f);
            y(beta_f_Nr) = R_far*cos(beta_f)*sin(alpha_f);
            z(beta_f_Nr) = R_far*sin(beta_f);
            if beta_f_Nr == length(beta_f_arr)
                text(x(beta_f_Nr),y(beta_f_Nr),1.15*z(beta_f_Nr), [num2str(alpha_f*180/pi) '^\circ'],'HorizontalAlignment','center','fontsize',6)
            end
        end
        plot3(x,y,z,'color', 'black','LineWidth',0.1)
    end

    beta_f_arr_meshlines = (-20:20:20)*pi/180;
    for beta_f_Nr = 1:length(beta_f_arr_meshlines)
        beta_f = beta_f_arr_meshlines(beta_f_Nr);
        x = zeros(size(alpha_f_arr));
        y = zeros(size(alpha_f_arr));
        z = zeros(size(alpha_f_arr));
        for alpha_f_Nr = 1:length(alpha_f_arr)
            alpha_f = alpha_f_arr(alpha_f_Nr);
            x(alpha_f_Nr) = R_far*cos(beta_f)*cos(alpha_f);
            y(alpha_f_Nr) = R_far*cos(beta_f)*sin(alpha_f);
            z(alpha_f_Nr) = R_far*sin(beta_f);
            if alpha_f_Nr == 1
                text(x(alpha_f_Nr),y(alpha_f_Nr)+0.05*x(alpha_f_Nr),z(alpha_f_Nr) - 0.02*R_far, [num2str(beta_f*180/pi) '^\circ'],'fontsize',6)
            end
        end
        plot3(x,y,z,'color', 'black', 'LineWidth',0.1)
    end
end
set(gca,'Color','none')
c = colorbar;
ax = gca;
axpos = ax.Position;
cpos = c.Position;
cpos(4) = 0.6*cpos(4);
cpos(3) = 0.3*cpos(3);
cpos(2) = 2.5*cpos(2);
cpos(1) = 0.85*cpos(1);
c.Position = cpos;
ax.Position = axpos;
if alpha_s_arr == 240*pi/180
    caxis([-40.738912799260937  24.031618135198975])
else
    caxis([-51.582570664320862  24.020551323538328])
end
export_fig(['../graphics/M3_SHBC_f' num2str(f_arr) '_' num2str(alpha_s_arr*180/pi) '_' num2str(max(beta_f_arr*180/pi))], '-png', '-transparent', '-r300')
