close all
clear all
% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot NURBS BeTSSi submarine
startMatlabPool
% setBCParameters
% M = 1;
% method = 'BEM';
% degreeElev = 0;
setBCParameters

P_inc = 1; % Amplitude of incident wave
rho_s = 7850; % Density of solid
rho_f = [1000, 1000]; % Density of fluids
c_f = [1500, 1500];  % Speed of sound in fluid domains
E = 210e9; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material

t = 0.01;   % Thickness

a = 7.0;
b = 3.5;
L = 42.0;

g2 = 6.5;
g3 = 6.5;
alpha = 18*pi/180;
beta = 240*pi/180;

c = 4.0;
s = 1.2;

b_ls = 2*s;
b_us = 2;
l_ls = 13;
l_us = 12.3;
h_s = 3.5;
delta_s = 0.2;
x_s = a-19;

l_lm = 2.6;
b_lm = 0.4;
l_um = 2.35;
b_um = 0.3;
h_m = 3.5;
delta_m = 0.25;
x_m = -51.9;


C_4 = c + b*cos(beta/2);
C_3 = b*sin(beta/2)-s;
C_2 = (2*C_4+C_3*tan(beta/2))/C_3^3;
C_1 = -(3*C_4+C_3*tan(beta/2))/C_3^2;

nPanels = 5;
dyPanels = C_3/nPanels;
P_panels = @(y) c + C_1*(y-s).^2 + C_2*(y-s).^3;

l_ld = 2.6;
b_ld = 2*(c-P_panels(s+dyPanels));
l_ud = 2.35;
b_ud = 0.22;
h_d = b-s;
delta_d = 0.25;
x_d = -4;


% t = b_lm/l_lm;
% t = b_ls/l_ls;
t = b_ld/l_ld;
f_lm = @(x) getNACA(sqrt(x),t,0);
aC = [0.2969; -0.1267; -0.3523; 0.2843; -0.1022]; % BeTSSi old
f_lsOld = @(x) 5*t*(aC(1)*sqrt(x)+aC(2)*x+aC(3)*x.^2+aC(4)*x.^3+aC(5)*x.^4);
integral(@(x)f_lm(x)-f_lsOld(x),0,1)/integral(@(x)f_lsOld(x),0,1)
% return
figure
scale = 4*80;
npts_xi = 1000;
xi = linspace(0,1,npts_xi).';
eta = linspace(0,1,2).';
xi_arr = copyVector(xi,2,1).';
eta_arr = copyVector(eta,npts_xi,2).';
% if 1   
%     h = g2*tan(alpha./2); % = g2./tan((pi-alpha)./2)
%     
%     x_c = -(L+g2+(b-h).*cot(alpha));
%     L_c = L+g2+(b-h).*cot(alpha);
%     x_f = -L-g2-g3/2;
%     P_C = [x_f;
%            (x_f+L_c)*tan(alpha)*cos(pi/2-30*pi/180);
%            (x_f+L_c)*tan(alpha)*sin(pi/2-30*pi/180)];
% 	
%     [g,dg,S,dSdxi,dSdeta] = getMainRudderFunctions(b,g2,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     P2 = @(xi) S(xi,g(xi));
% %     dP2 = @(xi) dSdxi(xi,g(xi)) + dSdeta(xi,g(xi)).*repmat(dg(xi),3,1);
% %     
% %     Eps = 1e-15;
% %     gTemp = @(xi) dot(dP2(xi),P2(xi)-P_C);
% %     xi1 = fzero(gTemp,0.184706062558860);
% %     return
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
% %     t = b_ls/l_ls;
%     xi1 = 0.184706062558860;
%     s2 = 1/8;
%     [~,~,arc] = NACAarcLength2(1,g,dg,dSdxi,dSdeta);
%     xi2 = invertNACA2_2(s2,l_lm,delta_m,x_m,g,dg);
%     totLength = arc(xi2);
%     s1 = NACAarcLength2(xi1,g,dg,dSdxi,dSdeta)*s2/totLength;
%     
%     S1 = getBeTSSiConePart(xi_arr,eta_arr,b,g2,g3,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m,xi1,s1,s2,'rudder',1);
%     r = getBeTSSiConePart(xi.',NaN,b,g2,g3,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m,xi1,s1,s2,'curve',1);
%     r2 = getBeTSSiConePart(NaN,eta.',b,g2,g3,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m,xi1,s1,s2,'curve2',1);
%     % 
%     npts2 = ceil(scale*s1);
%     npts3 = ceil(scale*(s2-s1));
%     npts4 = ceil(scale*(1-s2));
%     npts5 = ceil(scale/4);
%     xi = [linspace(0,s1,npts2).'; linspace2(s1,s2,npts3).'; linspace(s2,1,npts4).'];
% %     xi = linspace(0,s1,npts3).';
% %     xi = linspace2(s1,s2,npts3).';
% %     xi = linspace(s2,1,npts4).';
%     nptsxi = numel(xi);
%     eta = linspace(0,1,npts5).';
%     xi_arr = copyVector(xi,npts5,1).';
%     eta_arr = copyVector(eta,nptsxi,2).';
%     S2 = getBeTSSiConePart(xi_arr,eta_arr,b,g2,g3,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m,xi1,s1,s2,'cone',1);
% 
% 
%     h_semp = g2*tan(alpha/2); % = g2/tan((pi-alpha)/2)
%     theta = linspace(0,2*pi,npts_xi);
%     x = linspace(-L-g2-g3,-L-g2-g3/2,npts_xi);
%     [X,THETA] = meshgrid(x,theta);
%     R = (X+L+g2+cot(alpha)*(b-h_semp))/cot(alpha);
%     surf(X,R.*cos(THETA),R.*sin(THETA), 'EdgeColor','none','LineStyle','none', 'FaceColor', getColor(1),'FaceAlpha',0.4)
%     hold on
%     surf(reshape(S1(1,:),npts_xi,2),reshape(S1(2,:),npts_xi,2),reshape(S1(3,:),npts_xi,2), 'EdgeColor','none','LineStyle','none', 'FaceColor', getColor(1))
% 
% %     surf(reshape(S2(1,:),nptsxi,npts5),reshape(S2(2,:),nptsxi,npts5),reshape(S2(3,:),nptsxi,npts5), 'EdgeColor','none','LineStyle','none', 'FaceColor', getColor(1))
%     surf(reshape(S2(1,:),nptsxi,npts5),reshape(S2(2,:),nptsxi,npts5),reshape(S2(3,:),nptsxi,npts5), 'FaceColor', getColor(1))
%     plot3(r(1,:),r(2,:),r(3,:),'red')
%     plot3(r2(1,:),r2(2,:),r2(3,:),'red')
%     axis equal
%     view(-147,34)
%     view(-180,90)
%     camlight
%     ax = gca;               % get the current axis
%     ax.Clipping = 'off';    % turn clipping off
% 
%     return
% else
% %     
%     [g,dg,S,dSdxi,dSdeta] = getDepthRudderFunctions(b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,'upper');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     dx_rudder = 0.05;
% %     P_C = [x_d+dx_rudder; s; c];
% %     
% %     P2 = @(xi) S(xi,g(xi));
% %     dP2 = @(xi) dSdxi(xi,g(xi)) + dSdeta(xi,g(xi)).*repmat(dg(xi),3,1);
% %     maxIter = 100;
% %     Eps = 1e-15;
% %     gTemp = @(xi) dot(dP2(xi),P2(xi)-P_C);
% %     xi1 = fzero(gTemp,0.201905012203230);
% %     return
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     xi1 = 0.201905012203230;
%     s2 = 0.15;
%     [~,~,arc] = NACAarcLength2_depth(1,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,'rudderupper');
%     xi2 = invertNACA2_depth_2(s2,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,'rudderupper');
%     totLength = arc(xi2);
%     s1 = NACAarcLength2_depth(xi1,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,'rudderupper')*s2/totLength;
%     
%     S1u = getBeTSSiPanelPart(xi_arr,eta_arr,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2,'rudderupper');
%     S1l = getBeTSSiPanelPart(xi_arr,eta_arr,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2,'rudderlower');
%     hold on
%     surf(reshape(S1u(1,:),npts,npts),reshape(S1u(2,:),npts,npts),reshape(S1u(3,:),npts,npts), 'EdgeColor','none','LineStyle','none', 'FaceColor', getColor(1))
%     surf(reshape(S1l(1,:),npts,npts),reshape(S1l(2,:),npts,npts),reshape(S1l(3,:),npts,npts), 'EdgeColor','none','LineStyle','none', 'FaceColor', getColor(1))
% %     return
%     ru = getBeTSSiPanelPart(xi.',NaN,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2,'curveupper');
%     rl = getBeTSSiPanelPart(xi.',NaN,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2,'curvelower');
%     % 
%     npts2 = 5;
%     npts3 = 10;
%     npts4 = 40;
%     npts5 = 9;
%     xi = [linspace(0,s1,npts2).'; linspace2(s1,s2,npts3).'; linspace(s2,1,npts4).'];
% %     xi = linspace(0,s1,npts3).';
% %     xi = linspace2(s1,s2,npts3).';
% %     xi = linspace(s2,1,npts4).';
% %     xi = linspace(0.3,1,npts4).';
%     nptsxi = numel(xi);
%     eta = linspace(0,1,npts5).';
%     xi_arr = copyVector(xi,npts5,1).';
%     eta_arr = copyVector(eta,nptsxi,2).';
%     S2 = getBeTSSiPanelPart(xi_arr,eta_arr,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2,'lower');
%     S3 = getBeTSSiPanelPart(xi_arr,eta_arr,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2,'upper');
% 
% %     surf(reshape(S2(1,:),nptsxi,npts5),reshape(S2(2,:),nptsxi,npts5),reshape(S2(3,:),nptsxi,npts5), 'EdgeColor','none','LineStyle','none', 'FaceColor', getColor(1))
% %     surf(reshape(S3(1,:),nptsxi,npts5),reshape(S3(2,:),nptsxi,npts5),reshape(S3(3,:),nptsxi,npts5), 'EdgeColor','none','LineStyle','none', 'FaceColor', getColor(1))
%     surf(reshape(S3(1,:),nptsxi,npts5),reshape(S3(2,:),nptsxi,npts5),reshape(S3(3,:),nptsxi,npts5), 'FaceColor', getColor(1))
%     surf(reshape(S2(1,:),nptsxi,npts5),reshape(S2(2,:),nptsxi,npts5),reshape(S2(3,:),nptsxi,npts5), 'FaceColor', getColor(1))
%     plot3(ru(1,:),ru(2,:),ru(3,:),'red')
%     plot3(rl(1,:),rl(2,:),rl(3,:),'red')
%     % plot3(r2(1,:),r2(2,:),r2(3,:),'red')
%     axis equal
%     view(-147,34)
%     ax = gca;               % get the current axis
%     ax.Clipping = 'off';    % turn clipping off
%     % camlight
%     return
% end


% p = 8;
% % S = @(xi,eta) getBeTSSiPanelPart(xi,eta,a,b,l_ls,l_us,b_ls,b_us,h_s,delta_o,beta,s,c,xi_t,'rudderlower');
% r = @(xi) getBeTSSiPanelPart(xi.',zeros(size(xi.')),a,b,l_ls,l_us,b_ls,b_us,h_s,delta_o,beta,s,c,xi_t,'curveupper').';
% Xi = [zeros(1,p+1),xi_t*ones(1,p),sqrt(0.3)*ones(1,p),linspace2(sqrt(0.3),1,3),ones(1,p+1)];
% nurbs = leastSquares1D(Xi,p,r,3);
% sqrt(integral(@(xi) norm(evaluateNURBS(nurbs,xi).'-r(xi)).^2,0,1,'ArrayValued',true)/integral(@(xi)norm(r(xi)).^2,0,1,'ArrayValued',true))
% plotNURBSvec(nurbs,100,1,[0,0,0])
% maxE = -inf;
% for i = 1:numel(xi)
%     newE = max(norm(evaluateNURBS(nurbs,xi(i)).'-r(xi(i))));
%     if maxE < newE
%         maxE = newE;
%     end
% end
% sqrt(integral(@(xi) norm(evaluateNURBS(nurbs,xi).'-f(xi)).^2,0,1,'ArrayValued',true)/integral(@(xi)norm(f(xi)).^2,0,1,'ArrayValued',true))



nurbsCol2 = getBC_modData(a, b, L, g2, g3, alpha, beta, c, s);

% add sail
nurbsCol2(8) = addBeTSSiFoil(l_ls, b_ls, l_us, b_us, h_s, delta_s, [x_s,0,c], 0, true);

% add rudders
nurbsCol2(9) = addBeTSSiFoil(l_lm, b_lm, l_um, b_um, h_m, delta_m, [x_m,0,0], 0, true);
nurbsCol2(10) = addBeTSSiFoil(l_lm, b_lm, l_um, b_um, h_m, delta_m, [x_m,0,0], pi/2, true);
nurbsCol2(11) = addBeTSSiFoil(l_lm, b_lm, l_um, b_um, h_m, delta_m, [x_m,0,0], pi, true);
nurbsCol2(12) = addBeTSSiFoil(l_lm, b_lm, l_um, b_um, h_m, delta_m, [x_m,0,0], -pi/2, true);

% add depthrudders
nurbsCol2(13) = addBeTSSiFoil(l_ld, b_ld, l_ud, b_ud, h_d, delta_d, [x_d,s,c-b_ld/2], -pi/2, true);
nurbsCol2(14) = addBeTSSiFoil(l_ld, b_ld, l_ud, b_ud, h_d, delta_d, [x_d,-s,c-b_ld/2], pi/2, true);

parms = [a,b,L,g2,c,s,b_ls,b_us,l_ls,l_us,h_s,delta_s,x_s,l_lm,b_lm,l_um,b_um,h_m,delta_m,x_m,l_ld,b_ld,l_ud,b_ud,h_d,delta_d,...
         c+h_s,L+g2,L+g2+g3,x_s-l_ls,x_d-l_ld,x_m-l_lm];
nurbsCol2 = cleanNURBS(nurbsCol2,parms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    nurbsCol3 = getBCinteriorData();
    nurbsCol4 = nurbsCol2;
    nurbsCol4(8) = insertKnotsInNURBS(nurbsCol4(8),{sqrt(0.3)*ones(1,8) []});
    nurbsCol4(13) = insertKnotsInNURBS(nurbsCol4(13),{sqrt(0.3)*ones(1,8) []});
    nurbsCol4(14) = insertKnotsInNURBS(nurbsCol4(14),{sqrt(0.3)*ones(1,8) []});
    nurbsCol4(3) = insertKnotsInNURBS(nurbsCol4(3),{[] (L-12-l_ls*0.3)/L});
    nurbsCol4(5) = insertKnotsInNURBS(nurbsCol4(5),{[] (L-12-l_ls*0.3)/L});
    nurbsCol4(3) = insertKnotsInNURBS(nurbsCol4(3),{[] (L-4-l_ld*0.3)/L});
    nurbsCol4(5) = insertKnotsInNURBS(nurbsCol4(5),{[] (L-4-l_ld*0.3)/L});
    nurbsCol4 = cleanNURBS(nurbsCol4,parms);


%     figure
%     plotPart = 'none';
%     plotPart = 'nurbsCol2';
    plotPart = 'nurbsCol3';
%     plotPart = 'parts';
%     plotPart = 'CAD_transition';
    % plotPart = 'nurbsCol4';
%     plotPart = 'sail';
%     plotPart = 'depthrudder';
    lineWidth = 0.5;
%     lineWidth = 3;
    switch plotPart
        case 'parts'
            nurbsCol2 = explodeNURBS(nurbsCol2,'eta');
            resolution = [20 20];
%             resolution = [200 200];
            for i = 1:numel(nurbsCol2)
                if i == 24
                    plotNURBSvec(nurbsCol2{i},{'resolution',resolution,'color','blue','LineWidth',lineWidth});
                elseif i == 9 || i == 12 || i == 26
                    plotNURBSvec(nurbsCol2{i},{'resolution',resolution,'color','red','LineWidth',lineWidth});
                else
                    plotNURBSvec(nurbsCol2{i},{'resolution',resolution,'LineWidth',lineWidth});
                end
            end
            view(122,17) % standard
            camlight(150-122,70-17)
        case 'nurbsCol2'
            for i = 1:numel(nurbsCol2)
                plotNURBSvec(nurbsCol2{i},{'resolution',[200 200]});
            end
            view(18,10) % standard
        case 'nurbsCol3'
            steelColor = [67, 70, 75]/255;
            plotNURBSvec(nurbsCol2,{'resolution',[20 20], 'alphaValue',0.5, 'color',steelColor});
            plotNURBSvec(nurbsCol3(10:17),{'resolution',[10 10], 'alphaValue',0.5,  'color',steelColor});
            plotNURBSvec(nurbsCol3(setdiff(1:numel(nurbsCol3),10:17)),{'resolution',[10 10]});
            view(18,10) % standard
        case 'nurbsCol4'
            for i = 1:numel(nurbsCol4)
                plotNURBSvec(nurbsCol4{i},{'resolution',[10 10]});
            end
            view(18,10) % standard
        case 'CAD_transition'
            for i = 1:7
                if i == 6
                    plotNURBSvec(nurbsCol2{i},{'resolution',[200 200],'color','red'});
                else
                    plotNURBSvec(nurbsCol2{i},{'resolution',[200 200]});
                end
            end
            view(-60,20) % CAD_transition
            % export_fig('../graphics/BCA/CAD_transition', '-png', '-transparent', '-r300')
        case 'CAD_bow'
            for i = 1:7
                if i == 6
                    plotNURBSvec(nurbsCol2{i},{'resolution',[200 200],'color','red'});
                else
                    plotNURBSvec(nurbsCol2{i},{'resolution',[200 200]});
                end
            end
            view(120,20) % CAD_bow
        case 'tailSection'
            view(-140,30) % BeTSSi_BC_tailSection
            plotNURBSvec(nurbsCol2{7},{'resolution',[200 200]});
    %         export_fig('../graphics/BCA/BeTSSi_BC_tailSection', '-png', '-transparent', '-r300')
            plotControlPts(nurbsCol2{7}, [1 0 0], [1 0 0])
    %         export_fig('../graphics/BCA/BeTSSi_BC_tailSection_cp', '-png', '-transparent', '-r300')
    end
    axis equal
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off
    axis off
    % camlight
    camproj('perspective')

    if false
%         export_fig('../graphics/BCA/BCA_parts', '-png', '-transparent', '-r300')
        campos([24.421052461257915  24.531587345495094  21.785991629964464])
        camtarget([-16.950958582344327   0.383349541923005   5.973368116024933])
        camup([0, 0, 1])
        camva(6.305982631880568)
%         export_fig('../graphics/BCA/sailPart', '-png', '-transparent', '-r300', '-nocrop')

        campos([5.936172814767877   9.491372804675381   7.811160253655557])
        camtarget([-6.288869563455764   1.852318496730710   3.403897795199539])
        camup([0, 0, 1])
        camva(5.299051619974428)
%         export_fig('../graphics/BCA/drPart', '-png', '-transparent', '-r300', '-nocrop')

        view(131.2,28.2)
        campos([-38.799571606557109  12.783396618870956  12.888252320001310])
        camtarget([-52.838017930414537   0.493581355272471   2.883887889519321])
        camup([0, 0, 1])
        camva(6.540700628845736)
%         export_fig('../graphics/BCA/rPart', '-png', '-transparent', '-r300', '-nocrop')
    end
    if false
        campos([232.9067336994086   160.5331441878917   94.6177077941945])
        camtarget([-24.0000         0    2.0000])
        camup([0, 0, 1])
        camva(7.605458951342262)
%         export_fig('../graphics/BCA/BCA_error_p2', '-png', '-transparent', '-r300')
        campos([24.104013502569813  32.234272601422795  28.5])
        camtarget([-15.968494011405312   0.531647793061419   6.6])
        camup([0, 0, 1])
        camva(6.305982631880568)
%         export_fig('../graphics/BCA/BCA_error_p2_sail', '-png', '-transparent', '-r300', '-nocrop')

        campos([5.936172814767877   9.491372804675381   7.811160253655557])
        camtarget([-6.18   1.852318496730710   3.403897795199539])
        camup([0, 0, 1])
        camva(5.299051619974428)
%         export_fig('../graphics/BCA/BCA_error_p2_sail_drPart', '-png', '-transparent', '-r300', '-nocrop')

        view(131.2,28.2)
        campos([-38.799571606557109  12.783396618870956  12.888252320001310])
        camtarget([-52.838017930414537   0.493581355272471   2])
        camup([0, 0, 1])
        camva(6.540700628845736)
%         export_fig('../graphics/BCA/BCA_error_p2_rPart', '-png', '-transparent', '-r300', '-nocrop')
    end
    % return
    % camlight
    % view(0,0)
    % camlight(-30,0)
    % view(0,90)
    
    
    nurbsCol2 = explodeNURBS(nurbsCol2,'xi');
    nurbsCol2 = explodeNURBS(nurbsCol2,'eta');
    nurbsCol4 = explodeNURBS(nurbsCol4,'xi');
    nurbsCol4 = explodeNURBS(nurbsCol4,'eta');
%     for i = 1:numel(nurbsCol2)
%         printNURBSToFile(nurbsCol2{i}, ['~/Onedrive/work/rhinoceros/BCdata/BeTSSi_mod_' num2str(i) '.txt'])
%     end
%     for i = 1:numel(nurbsCol4)
%         printNURBSToFile(nurbsCol4{i}, ['~/Onedrive/work/rhinoceros/BCdata/BeTSSi_mod2_' num2str(i) '.txt'])
%     end
%     for i = 1:numel(nurbsCol3)
%         printNURBSToFile(nurbsCol3{i}, ['~/Onedrive/work/rhinoceros/BCdata/BeTSSi_interior' num2str(i) '.txt'])
%     end
%     
%     for p = [2,3,4]
%         load(['NURBSgeometries/BCAdata/BeTSSi_BCA_p' num2str(p)])
%         for i = 1:numel(nurbs)
%             printNURBSToFile(nurbs{i}, ['~/Onedrive/work/rhinoceros/BCdata/BeTSSi_mod_p' num2str(p) '_' num2str(i) '.txt'])
%         end
%     end

    if strcmp(plotPart, 'nurbsCol2') || strcmp(plotPart, 'parts')
        [g,dg,S,dSdxi,dSdeta] = getMainRudderFunctions(b,g2,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m);

        xi1 = 0.184706062558860;
        s2 = 1/8;
        [~,~,arc] = NACAarcLength2(1,g,dg,dSdxi,dSdeta);
        xi2 = invertNACA2_2(s2,l_lm,delta_m,x_m,g,dg);
        totLength = arc(xi2);
        s1 = NACAarcLength2(xi1,g,dg,dSdxi,dSdeta)*s2/totLength;

        npts_xi = 1000;
        xi = linspace(0,1,npts_xi).';
        r = getBeTSSiConePart(xi.',NaN,b,g2,g3,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m,xi1,s1,s2,'curve',1);
        r2 = r;
        r(2,:) = -r(2,:);
        r = [r,r2(:,end:-1:1)];
        for theta = [0,pi/2,pi,3*pi/2]
            R = rotationMatrix(theta, 'Xaxis');
            r = R*r;
            plot3(r(1,:),r(2,:),r(3,:),'color','black','LineWidth',lineWidth)
        end

        [g,dg,S,dSdxi,dSdeta] = getDepthRudderFunctions(b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,'upper');

        xi1 = 0.201905012203230;
        s2 = 0.15;
        [~,~,arc] = NACAarcLength2_depth(1,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,'rudderupper');
        xi2 = invertNACA2_depth_2(s2,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,'rudderupper');
        totLength = arc(xi2);
        s1 = NACAarcLength2_depth(xi1,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,'rudderupper')*s2/totLength;

        ru = getBeTSSiPanelPart(xi.',NaN,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2,'curveupper');
        rl = getBeTSSiPanelPart(xi.',NaN,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2,'curvelower');
        r = [ru,rl(:,end:-1:1)];
        plot3(r(1,:),r(2,:),r(3,:),'color','black','LineWidth',lineWidth)
        r(2,:) = -r(2,:);
        plot3(r(1,:),r(2,:),r(3,:),'color','black','LineWidth',lineWidth)
    %     export_fig('../../graphics/BCA/CAD', '-png', '-transparent', '-r300')
%         export_fig('../../graphics/BCA/BCA_topView', '-png', '-transparent', '-r300')
%         savefig('../results/geometries/BeTSSi')
    end
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% varCol.colorFun = @(x) log10(minDistNURBS(x,nurbsCol2));
% nurbsCol2 = explodeNURBS(nurbsCol2{7},'eta');
nurbsCol2 = explodeNURBS(nurbsCol2{7},'xi');

indices = ones(323,1);
indices(1:16) = 1; % back disk
indices([17:32,33:37,  53:57,58:63,64:66,  76:78,79:80,81:83,  93:95,96:97,98:100,  110:112,113:118,  119:138]) = 2; % cone parts
indices(38:42) = 44; % rudder top
indices(43:47) = 45; % rudder top
indices(48:52) = 46; % rudder top
indices(67:69) = 47; % rudder righ_s
indices(70:72) = 48; % rudder righ_s
indices(73:75) = 49; % rudder righ_s
indices(84:86) = 50; % rudder bottom
indices(87:89) = 51; % rudder bottom
indices(90:92) = 52; % rudder bottom
indices(101:103) = 53; % rudder left
indices(104:106) = 54; % rudder left
indices(107:109) = 55; % rudder left

indices(139:150) = 4:15; % transition part
indices(151:158) = 3;   % lower transition part
% 
indices(159:167) = 17; % 1st left side panel
indices(168:176) = 18; % 2st left side panel
indices(177:185) = 19; % 3rd left side panel
indices(186:2:195) = 20; % 4th left side panel
indices(187:2:195) = 21; % 5th left side panel
indices(196:199) = 20; % 4th left side panel rudder
indices(200:203) = 56; % left depth rudder
indices(204:207) = 57; % left depth rudder
indices(208:211) = 58; % left depth rudder
indices(212:215) = 21;
indices(216) = 20;
indices(217) = 21;
indices([218:219,220:223, 236:239,240:249]) = 22; % deck
indices(224:227) = 41; % sail
indices(228:231) = 42; % sail
indices(232:235) = 43; % sail

indices(250:2:259) = 23;
indices(251:2:259) = 24;
indices(260:263) = 23;
indices(264:267) = 59; % righ_s depth rudder
indices(268:271) = 60; % righ_s depth rudder
indices(272:275) = 61; % righ_s depth rudder
indices(276:279) = 24;
indices(280) = 23;
indices(281) = 24;
indices(282:290) = 25;
indices(291:299) = 26;
indices(300:308) = 27;
indices(309:380) = 16; %lower cylindrical part
indices(381:392) = [29:34,34:39];
indices(393:400) = 28; % lower part of bow
% 
% conn = NaN(323,4);
% idx = 1:16;
% conn(idx,1) = idx(2:end)
% return
pArr = 2:9;
% pArr = 2;
% pArr = 15;
maxCp = zeros(size(pArr));
for ip = 1:numel(pArr)
    p = pArr(ip);
    if 0
        nurbsColApprox = getBCAData(p);
        nurbsColApprox = cleanNURBS(nurbsColApprox,parms);
        
        nurbs = nurbsColApprox;
%         save(['NURBSgeometries/BCAdata/BeTSSi_BCA_p' num2str(p)],'nurbs')
    else
        load(['NURBSgeometries/BCAdata/BeTSSi_BCA_p' num2str(p)])
        nurbsColApprox = nurbs;
    end
    tic
%     figure(36)
%     resolution = [80 80];
%     maxC = -Inf(numel(nurbsColApprox),1);
    I = 0;
    A = 0;
    for i = 1:numel(nurbsColApprox) %[1,2,4,6,8,10,11,13] % numel(nurbsColApprox) = 202
        temp_i = indices(i) - 3;
        if temp_i > 12 || temp_i < 1
            temp_i = 1;
        end
        fprintf([num2str(i), ', '])
        if ~mod(i,10)
            fprintf('\n')
        end
        f = @(v) uTest_BC(v,a,b,c,L,g2,g3,alpha,beta,s,indices(i),nurbsCol2{temp_i});
%         colorFun = @(v) log10(f(v));
%         [~,maxC(i)] = plotNURBSvec(nurbsColApprox{i},{'resolution',resolution, 'colorFun',colorFun,'elementBasedSamples',true, ...
%                                                     'samplingDistance',0.01,'LineWidth',3});
        [Itemp,Atemp] = computeGeometricL2(nurbsColApprox{i},{'f',f,'extraGP',10});
        I = I + Itemp;
        A = A + Atemp;
%         maxC(i)
%         if maxC(i) > -3
%             keyboard
%         end
    end
    maxCp(ip) = sqrt(I)/A;
    maxCp(ip)
%     axis equal
%     caxis([-15,-2])
% %     view(18,10)
%     view(150,45)
%     colorbar off
%     % view(180,90)
%     % camlight
%     ax = gca;               % get the current axis
%     ax.Clipping = 'off';    % turn clipping off
%     axis off
% %     plotControlPts(nurbsColApprox(33:47)) 
%     drawnow
%     if 1
% %         savefig(['results/BCA/geometryErrors/BeTSSi_BCA_p' num2str(p)]) % Max Error = 10^(-2.53) = 0.00295 
%         savefig(['results/BCA/geometryErrors/BeTSSi_BCA_p' num2str(p) '_linewidth']) % Max Error = 10^(-2.53) = 0.00295  
%         close all      
%     end
end
semilogy(maxCp)
% save('results/BCA/convergence','maxCp')
% printResultsToFile2('results/BCA/convergence', 10.^pArr.', maxCp.')
L_gamma = a+L+g2+g3;
printResultsToFile2('results/BCA/convergence2', pArr.', 100*maxCp.'/L_gamma)
return
for p = 2:15
    uiopen(['results/BCA/geometryErrors/BeTSSi_BCA_p' num2str(p) '.fig'],1)
    caxis([-15,-2])
    campos([232.9067336994086   160.5331441878917   94.6177077941945])
    camtarget([-24.0000         0    2.0000])
    camup([0, 0, 1])
    camva(7.605458951342262)
    savefig(['results/BCA/geometryErrors/BeTSSi_BCA_p' num2str(p)]) % Max Error = 10^(-2.53) = 0.00295  
    close all      
end
npts = 1000;
xi = linspace(0,1,npts);
eta = linspace(0,1,npts);
xi_arr = copyVector(xi,npts,1);
eta_arr = copyVector(eta,npts,2);
XX = X([xi_arr.'; eta_arr.']);
surf(reshape(XX(1,:).',npts,npts),reshape(XX(2,:).',npts,npts),reshape(XX(3,:).',npts,npts), 'EdgeColor','none','LineStyle','none', 'FaceColor', getColor(1),'FaceAlpha',0.5)

[XI,ETA] = meshgrid(xi,eta);


f_u = @(xi) getNACA(xi,b_ls/l_ls,0);
dfxi_u = @(xi) getNACA(xi,b_ls/l_ls,1);
dfdxi2_u = @(xi) getNACA(xi,b_ls/l_ls,2);
f_o = @(xi) getNACA(xi,b_us/l_us,0);
dfxi_o = @(xi) getNACA(xi,b_us/l_us,1);
dfdxi2_o = @(xi) getNACA(xi,b_us/l_us,2);
surf(l_ls*XI.^2+ETA.*(delta_s+(l_us-l_ls)*XI.^2), l_ls*f_u(XI) + ETA.*(l_us*f_o(XI)-l_ls*f_u(XI)),ETA*h_s, 'EdgeColor','none','LineStyle','none', 'FaceColor', getColor(1))

for i = 1:numel(nurbsColApprox)
    degree = nurbsColApprox{i}.degree;
    if any(degree < 2)
        nurbsColApprox{i} = elevateNURBSdegree(nurbsColApprox{i},[2,2]-degree);
    end
end
for i = 1:numel(nurbsColApprox)
    P = nurbsColApprox{i}.coeffs(1:3,2,2);
    text(P(1),P(2),P(3),num2str(i),'Color','red')
end

%%%%%%%%%%%%%%%%%%%%%%%
figure(42)
for i = 1:numel(nurbsCol2)
    plotNURBSvec(nurbsCol2{i},[10 10], 1, getColor(1));
end
for i = 1:numel(nurbsCol2)
    degree = nurbsCol2{i}.degree;
    if any(degree < 2)
        elevDegree = [2,2]-degree;
        elevDegree(elevDegree < 0) = 0;
        nurbsCol2{i} = elevateNURBSdegree(nurbsCol2{i},elevDegree);
    end
end
for i = 1:numel(nurbsCol2)
    P = nurbsCol2{i}.coeffs(1:3,2,2);
    text(P(1),P(2),P(3),num2str(i))
end
view(18,10)
axis equal

%%%%%%%%%%%%%%%%%%%%%%%

% plotNURBSvec(nurbsCol2{1},[100 100], 1, getColor(1));
% hold on
% for i = 1:numel(nurbsCol2E)
%     varCol.colorFun = @(x) norm(x);
%     plotNURBSvec(nurbsCol2E{i},[10 10], 1, getColor(1));
% end
% axis equal
% view(122,16)
% % view(180,90)
% camlight
% axis off
% plotControlPts(nurbsCol{3})
% % [solid, hyp] = getBCData(parms);
% fluid = extractOuterSurface(solid);
% 
% 
% Xi = fluid.knots{1};
% h = g2*tan(alpha/2); % = g2/tan((pi-alpha)/2)
% x2 = g3*tan(alpha);
% 
% plotNURBSvec(fluid,[100 100], 1, getColor(1), 1);
% 
% % varCol.colorFun = @(v,xi,eta) uTest_BC(v,xi,eta,fluid.knots{1},fluid.knots{2},a,b,l,g2,g3,alpha);
% % plotNURBSvec(fluid,[100 100], 1, getColor(1), 1, NaN, varCol);
% 
% axis equal
% axis off
% set(gca, 'Color', 'none');
% % view(-70,30)
% view(120,10)
% drawnow
% % colorbar
% camlight
% % % camproj('perspective')
% % hold on
% % % plotControlPts(solid)
% % return
% 
% % 
% % % export_fig('../graphics/BeTSSi_BC_mod', '-png', '-transparent', '-r300')
% % % export_fig('../graphics/BeTSSi_BC_stripped', '-png', '-transparent', '-r300')
% % 
% 
% return
% 
% % export_fig('../graphics/BeTSSi_BC', '-png', '-transparent', '-r300')
% 
% % [nurbs, nurbs2, nurbs3, nurbs4, nurbs5] = getBCDataParts(a, b, L, g2, g3, alpha, beta, s, c);
% % varCol.colorFun = @(v) norm(v(2:3));
% % plotNURBSvec(nurbs,[100 100], 1, getColor(1), 1);
% % hold on
% % plotNURBSvec(nurbs2,[100 100], 1, [1 0 0], 1);
% % axis equal
% % axis off
% % set(gca, 'Color', 'none');
% % % view(-70,30)
% % view(20,36)
% % drawnow
% % camlight
% % camproj('perspective')
% % hold on
% % plotControlPts(nurbs)
% % figure(43)
% % plotNURBSvec(nurbs,[100 100], 1, getColor(1), 1);
% % hold on
% % plotNURBSvec(nurbs2,[100 100], 1, getColor(1), 1);
% % plotNURBSvec(nurbs3,[0 0], 1, [1 0 0], 1);
% % plotNURBSvec(nurbs4,[0 0], 1, [1 0 0], 1);
% % plotNURBSvec(nurbs5,[0 0], 1, [0 0 1], 1);
% % axis equal
% % axis off
% % set(gca, 'Color', 'none');
% % % view(-70,30)
% % view(55,34)
% % drawnow
% % camlight
% % camproj('perspective')
% % hold on
% % x0 = 0;
% % y0 = 0;
% % z0 = 0;
% % ss = 2;
% % quiver3(0,0,0,0,0,ss,'color','black','MaxHeadSize',0.3)
% % quiver3(0,0,0,0,ss,0,'color','black','MaxHeadSize',0.3)
% % quiver3(0,0,0,ss,0,0,'color','black','MaxHeadSize',0.3)
% % text([x0+ss, x0, x0-0.05*ss], [y0, y0+ss, y0-0.05*ss], [z0, z0, z0+1.1*ss], ['$$x$$';'$$y$$';'$$z$$'],'interpreter','latex');
% % 
% % figure(44)
% % [nurbs6, nurbs7,  nurbs8,  nurbs9,  nurbs10] = getBCDataParts2(a, b, L, g2, g3, alpha, beta, s, c);
% % if 1
% %     plotNURBSvec(nurbs,[200 400], 1, getColor(1), 1);
% %     plotNURBSvec(nurbs3,[0 0], 1, getColor(1), 1);
% %     plotNURBSvec(nurbs4,[0 0], 1, getColor(1), 1);
% %     plotNURBSvec(nurbs5,[0 0], 1, getColor(1), 1);
% %     plotNURBSvec(nurbs2,[100 100], 1, getColor(1), 1);
% %     % plotNURBSvec(nurbs8,[100 100], 1, getColor(1), 1);
% %     % plotNURBSvec(nurbs9,[100 100], 1, getColor(1), 1);
% %     % plotNURBSvec(nurbs10,[100 100], 1, getColor(1), 1);
% % 
% % 
% %     if 0
% %         plotNURBSvec(nurbs6,[100 100], 1, getColor(1), 1);
% %         plotNURBSvec(nurbs7,[100 100], 1, [1 0 0], 1);
% % 
% %         view(120,20)
% %     %     export_fig('../graphics/BeTSSi_BCpart4', '-png', '-transparent', '-r300')
% %     elseif 0
% %         view(300,20)
% %         plotNURBSvec(nurbs6,[100 100], 1, [1 0 0], 1);
% %         plotNURBSvec(nurbs7,[100 100], 1, getColor(1), 1);
% %     %     export_fig('../graphics/BeTSSi_BCpart3', '-png', '-transparent', '-r300')
% %     else
% %         view(18,10)
% %         plotNURBSvec(nurbs6,[200 400], 1, getColor(1), 1);
% %         plotNURBSvec(nurbs7,[200 400], 1, getColor(1), 1);
% %     %     export_fig('../graphics/BC/BeTSSi_BC_stripped', '-png', '-transparent', '-r600')
% % 
% %     end
% % else
% %     plotNURBSvec(nurbs6,[200 200], 1, getColor(1), 1);
% % %     view(30,30)
% %     view(-140,30)
% %     %     export_fig('../graphics/BeTSSi_BC_tailSection', '-png', '-transparent', '-r300')
% % %     plotControlPts(nurbs6, [1 0 0], [1 0 0])
% %     %     export_fig('../graphics/BeTSSi_BC_tailSection_cp', '-png', '-transparent', '-r300')
% % end
% % axis equal
% % axis off
% % set(gca, 'Color', 'none');
% % drawnow
% % camproj('perspective')
% % hold on
% %     
% % camlight
%  
% % figure(45)
% % [nurbs6, nurbs7,  nurbs8,  nurbs9,  nurbs10] = getBCDataParts2(a, b, l, g2, g3, alpha, beta, s, c);
% % plotNURBSvec(nurbs,[200 400], 1, getColor(1), 1);
% % axis equal
% % axis off
% % set(gca, 'Color', 'none');
% % % view(-70,30)
% % view(55,34)
% % drawnow
% % camproj('perspective')
% % hold on
% % plotNURBSvec(nurbs3,[0 0], 1, getColor(1), 1);
% % plotNURBSvec(nurbs4,[0 0], 1, getColor(1), 1);
% % plotNURBSvec(nurbs5,[0 0], 1, getColor(1), 1);
% % plotNURBSvec(nurbs2,[100 100], 1, getColor(1), 1);
% % plotNURBSvec(nurbs6,[100 100], 1, getColor(1), 1);
% % plotNURBSvec(nurbs7,[100 100], 1, getColor(1), 1);
% % 
% % view(18,10)
% % camlight
% % 
% % 
% % hold on
% % % add sail
% % addSailToBeTSSi_exact(a0, a1, a2, a3, a4, l_ls, b_ls, l_us, b_us, h_s, delta_s, [a-19.0,0,c], 0);
% % l_ls = 2.6;
% % b_ls = 0.4;
% % l_us = 2.35;
% % b_us = 0.3;
% % delta_s = 0.25;
% % % add rudders
% % addSailToBeTSSi_exact(a0, a1, a2, a3, a4, l_ls, b_ls, l_us, b_us, h_s, delta_s, [a-58.9,0,0], pi/2);
% % addSailToBeTSSi_exact(a0, a1, a2, a3, a4, l_ls, b_ls, l_us, b_us, h_s, delta_s, [a-58.9,0,0], 0);
% % addSailToBeTSSi_exact(a0, a1, a2, a3, a4, l_ls, b_ls, l_us, b_us, h_s, delta_s, [a-58.9,0,0], -pi/2);
% % addSailToBeTSSi_exact(a0, a1, a2, a3, a4, l_ls, b_ls, l_us, b_us, h_s, delta_s, [a-58.9,0,0], pi);
% % 
% % b_ls = 0.265;
% % b_us = 0.22;
% % h_s = b-s;
% % % add depthrudders
% % addSailToBeTSSi_exact(a0, a1, a2, a3, a4, l_ls, b_ls, l_us, b_us, h_s, delta_s, [a-11.0,s,c-0.265], -pi/2);
% % addSailToBeTSSi_exact(a0, a1, a2, a3, a4, l_ls, b_ls, l_us, b_us, h_s, delta_s, [a-11.0,-s,c-0.265], pi/2);
% 
% % % %     export_fig('../graphics/BeTSSi_BC', '-png', '-transparent', '-r600')
