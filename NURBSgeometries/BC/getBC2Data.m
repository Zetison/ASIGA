function [nurbsCol, hyp] = getBC2Data(parms)

a = parms.a;
b = parms.b;
L = parms.L;
g2 = parms.g2;
g3 = parms.g3;
alpha = parms.alpha;
c = parms.c;
s = parms.s;
t = parms.t;
l_ls = parms.l_ls;
l_us = parms.l_us;
b_ls = parms.b_ls;
b_us = parms.b_us;
h_s = parms.h_s;
delta_s = parms.delta_s;
[nurbsCol, hyp] = getBC2Data_0(a, b, L, g2, g3, alpha, c, s, l_ls, b_ls, l_us, b_us, h_s, delta_s);
% nurbsCol = getBC2Data_0(a-t, b-t, L, g2-t/2, g3-t/2, alpha, c-t, s-t/2, l_ls, b_ls, l_us, b_us, h_s, delta_s);
% Zeta = [0 0 1 1];
% coeffs = nurbs0.coeffs;
% coeffs(:,:,:,2) = nurbs1.coeffs;

% nurbs = createNURBSobject(coeffs,{nurbs1.knots{1}, nurbs1.knots{2}, Zeta});

function [nurbsCol, hyp] = getBC2Data_0(a, b, L, g2, g3, alpha, c, s, l_ls, b_ls, l_us, b_us, h_s, delta_s)

t = b_ls/l_ls;
N1 = 1;
N2 = 6;
p = 2;
xi_T = 0.168967783470083;
[cntrlPts, Xi] = getNACAapprox(t,p,N1,N2,xi_T);
cntrlPts = flipud(cntrlPts); 
cntrlPts(:,1) = 1-cntrlPts(:,1); 
h = g2*tan(alpha/2); % = g2/tan((pi-alpha)/2)

x1 = tan(alpha/2)*(g2/tan(alpha) + h);
x2 = g3*tan(alpha);

lowerCpts = [  -L-g2-g3,    0,  0,          1; % 0
               -L-g2-g3,    0, (-b+h+x2)/2, 1;
               -L-g2-g3,    0, -b+h+x2,     1; % 1
               -L-g2-g3/2,  0, -b+h+x2/2,   1;
               -L-g2,       0, -b+h,        1; % 2
               -L-x1,       0, -b,          cos(alpha/2);
               -L,          0, -b,          1; % 3
               -L/2,        0, -b,          1;
                0,          0, -b,          1; % 4
                a,          0, -b,          1/sqrt(2);
                a,          0,  0,          1];% 5
            
righ_sCpts1 = lowerCpts;
righ_sCpts1(:,2) = abs(lowerCpts(:,3))*tan(pi/4);
righ_sCpts1(:,4) = righ_sCpts1(:,4)*cos(pi/4);



R_x = rotationMatrix(pi/2, 'Xaxis');

righ_sCpts2 = lowerCpts;
righ_sCpts2(:,1:3) = (R_x*righ_sCpts2(:,1:3).').';

upperCpts = lowerCpts;
upperCpts(:,3) = abs(upperCpts(:,3));

righ_sCpts3 = upperCpts;
righ_sCpts3(:,2) = abs(upperCpts(:,3))*tan(pi/2-pi/4);
righ_sCpts3(:,4) = righ_sCpts3(:,4)*cos(pi/2-pi/4);

coeffs = zeros(4,5,11);
coeffs(:,1,:) = lowerCpts.';
coeffs(:,2,:) = righ_sCpts1.';
coeffs(:,3,:) = righ_sCpts2.';
coeffs(:,4,:) = righ_sCpts3.';
coeffs(:,5,:) = upperCpts.';


Xi = [0 0 0 0.5 0.5 1 1 1];
Eta = [0 0 0 1 1 2 2 3 3 4 4 5 5 5]/5;

nurbs = createNURBSobject(coeffs,{Xi, Eta});
xi_r = 0.80;
nurbs = insertKnotsInNURBS(nurbs,{[xi_r xi_r 0.875] 0.5});
nurbs.coeffs(3,end-2:end,6:end-1) = 4;

scale = (0.875 - xi_r)/0.5;
nurbs.coeffs(1:3,5,7:end) = (1-scale)*nurbs.coeffs(1:3,6,7:end)+scale*nurbs.coeffs(1:3,4,7:end);
coeffs = nurbs.coeffs;
for i = 4:size(coeffs,2)    
    v1 = coeffs(1:3,i,4);
    v2 = coeffs(1:3,i,5);
    
    x_g4 = coeffs(1,1,6);
    
    coeffs(1:3,i,6) = v1 + (x_g4-v1(1))/(v2(1)-v1(1))*(v2-v1);
end
nurbs.coeffs = coeffs;

%% create hole for sail
eta_s = 0.681; 
eta_e = 0.75;
newKnots = @(eta_s,eta_e) [eta_s eta_s (eta_s+cntrlPts(2:end-2,1).'*(eta_e-eta_s)) eta_e eta_e];
nurbs = insertKnotsInNURBS(nurbs,{[] newKnots(eta_s,eta_e)});
nurbs.coeffs(1,end,10:19) = -19+a-l_ls+l_ls*cntrlPts(:,1);
nurbs.coeffs(2,end,10:19) = l_ls*cntrlPts(:,2);
nurbs.coeffs(3,end,10:19) = c;
nurbs.coeffs(4,end,10:19) = 1;
nurbs.coeffs(1:3,end-1,7:end) = (nurbs.coeffs(1:3,end,7:end)+nurbs.coeffs(1:3,end-2,7:end))/2;

P1 = nurbs.coeffs(1:3,6,18);
P2 = nurbs.coeffs(1:3,6,19);
P3 = nurbs.coeffs(1:3,7,19);
nurbs.coeffs(1:3,7,18) = P3 + P1-P2;

P1 = nurbs.coeffs(1:3,end,11);
P2 = nurbs.coeffs(1:3,end-2,11);
P3 = nurbs.coeffs(1:3,end-1,10);
nurbs.coeffs(1:3,end-1,11) = P1 + (P3(2)-P1(2))/(P2(2)-P1(2))*(P2-P1);

%% sail
coeffs = ones(4,4,10);
Xi = [0,0,1,2,3,3]/3;
Eta = (nurbs.knots{2}(11:end-5)-eta_s)/(eta_e-eta_s);
Eta = [0,Eta,1];
coeffs(1:3,1,:) = nurbs.coeffs(1:3,end,10:19);
coeffs(1,2,:) = -delta_s-19+a-l_us+l_us*cntrlPts(:,1);
coeffs(2,2,:) = l_us*cntrlPts(:,2);
coeffs(3,2,:) = c+h_s;

coeffs(:,3:4,:) = coeffs(:,2:-1:1,:);
coeffs(2,3:4,:) = -coeffs(2,3:4,:);

nurbs_sail = createNURBSobject(coeffs,{Xi, Eta});
nurbs_sail = elevateNURBSdegree(nurbs_sail,[1,0]);
temp = nurbs_sail.coeffs(1,4,end-4:end-1);
nurbs_sail.coeffs(1,4,end-4:end-1) = temp(:) - l_ls*[0.005; 0.01; 0.023; 0.031];

%% create hole for depth rudders
eta_s = 0.765;
eta_e = 0.785;
l_ls = 2.6;
l_us = 2.35;
delta_s = l_ls-l_us;
slope = 2.5;
nurbs = insertKnotsInNURBS(nurbs,{[xi_r*0.95, xi_r] newKnots(eta_s,eta_e)});
nurbs.coeffs(1,end-3,21:30) = -11+a-l_ls+l_ls*cntrlPts(:,1);
nurbs.coeffs(2,end-3,21:30) = nurbs.coeffs(2,end-3,21) - cntrlPts(:,2)*slope;
nurbs.coeffs(3,end-3,21:30) = nurbs.coeffs(3,end-3,20)+l_ls*cntrlPts(:,2);
nurbs.coeffs(4,end-3,21:30) = 1;

nurbs.coeffs(1,end-4,21:30) = -11+a-l_ls+l_ls*cntrlPts(:,1);
nurbs.coeffs(2,end-4,21:30) = nurbs.coeffs(2,end-4,21) + cntrlPts(:,2)*slope;
nurbs.coeffs(3,end-4,21:30) = nurbs.coeffs(3,end-4,20)-l_ls*cntrlPts(:,2);
nurbs.coeffs(4,end-4,21:30) = 1;

nurbs.coeffs(1,end,21:30) = nurbs.coeffs(1,1,21:30);
nurbs.coeffs(1,end-1,21:30) = nurbs.coeffs(1,1,21:30);

%% depth rudder
coeffs = ones(4,4,10);
Xi = [0,0,1,2,3,3]/3;
Eta = (nurbs.knots{2}(22:end-5)-eta_s)/(eta_e-eta_s);
Eta = [0,Eta,1];
coeffs(1:3,1,:) = nurbs.coeffs(1:3,end-4,21:30);
coeffs(1,2,:) = -delta_s-11+a-l_us+l_us*cntrlPts(:,1);
coeffs(3,2,:) = nurbs.coeffs(3,end-4,21)-l_us*cntrlPts(:,2);
coeffs(1,3,:) = -delta_s-11+a-l_us+l_us*cntrlPts(:,1);
coeffs(3,3,:) = nurbs.coeffs(3,end-4,21)+l_us*cntrlPts(:,2);
coeffs(2,2:3,:) = b;
coeffs(1:3,4,:) = nurbs.coeffs(1:3,end-3,21:30);
nurbs_rudder = createNURBSobject(coeffs,{Xi, Eta});
nurbs_rudder = elevateNURBSdegree(nurbs_rudder,[1,0]);
temp = nurbs_rudder.coeffs(1,4,end-4:end-1);
nurbs_rudder.coeffs(1,4,end-4:end-1) = temp(:) - l_us*[0.005; 0.01; 0.023; 0.031];

nurbs_rudder2 = nurbs_rudder;
nurbs_rudder2.coeffs(2,:,:) = -nurbs_rudder2.coeffs(2,:,:);
% nurbs_rudder = insertKnotsInNURBS(nurbs_rudder,{[linspace2(0,1/3,3),linspace2(1/3,2/3,1),linspace2(2/3,1,3)] []});


%% create hole for back rudders
eta_s = 0.215;
eta_e = 0.31;
h_s = 3.5;
nurbs = insertKnotsInNURBS(nurbs,{[0.125, 0.375, 0.5, 0.625] newKnots(eta_s,eta_e)});
nurbs.coeffs(1,6,5) = -58.9+a-l_ls;
nurbs.coeffs(2,6,5) = abs(-b+h+x2) + (L+g2+g3+a-58.9-l_ls)*tan(alpha);
nurbs.coeffs(1,6,6) = -58.9+a-l_ls + l_ls*cntrlPts(2,1);
coeffs = nurbs.coeffs;
for i = 1:size(coeffs,2)
    v1 = coeffs(1:3,i,3);
    v2 = coeffs(1:3,i,4);
    
    x_g4 = coeffs(1,6,5);
        
    coeffs(1:3,i,5) = v1 + (x_g4-v1(1))/(v2(1)-v1(1))*(v2-v1);
    coeffs(1:3,i,4) = (coeffs(1:3,i,3)+coeffs(1:3,i,5))/2;
    
    v2 = coeffs(1:3,i,4);
    v3 = coeffs(1:3,i,5);
    x_g5 = coeffs(1,6,6);
    coeffs(1:3,i,6) = v2 + (x_g5-v2(1))/(v3(1)-v2(1))*(v3-v2);
end
nurbs.coeffs = coeffs;


nurbs.coeffs(1,6,6:14) = -58.9+a-l_ls+l_ls*cntrlPts(2:end,1);
nurbs.coeffs(2,6,14) = abs(-b+h+x2) + (L+g2+g3+a-58.9)*tan(alpha);
y_1 = nurbs.coeffs(2,6,5);
y_2 = nurbs.coeffs(2,6,14);
nurbs.coeffs(2,6,6:13) = y_1 + (y_2-y_1)*cntrlPts(2:end-1,1);

nurbs.coeffs(3,6,6:14) = 0+l_ls*cntrlPts(2:end,2);
nurbs.coeffs(4,6,5:14) = 1;
temp = nurbs.coeffs(3,7,6:13);
nurbs.coeffs(3,7,6:13) = temp(:) + 0.4*l_ls*cntrlPts(2:end-1,2);

nurbs.coeffs(:,5,5:14) = nurbs.coeffs(:,6,5:14);
nurbs.coeffs(3,5,5:14) = -nurbs.coeffs(3,6,5:14);
nurbs.coeffs(:,4,5:14) = nurbs.coeffs(:,7,5:14);
nurbs.coeffs(3,4,5:14) = -nurbs.coeffs(3,7,5:14);

%% Extract temp
% Xi = [0, (nurbs.knots{1}(7:13)-0.5)/(xi_r-0.5)];
% Eta = [nurbs.knots{2}(1:16)/eta_e, 1];
% nurbs_temp = createNURBSobject(nurbs.coeffs(:,6:10,1:14),{Xi, Eta});
% 
% close all
% varCol.colorFun = @(v,xi,eta) log10(uTest_BC(v,xi,eta,nurbs_temp.knots{1},nurbs_temp.knots{2},a,b,L,g2,g3,alpha,[2,numel(unique(Eta)),numel(unique(Eta))]));
% plotNURBS(nurbs_temp,[100 100], 1, 1.5*[44 77 32]/255, 1, NaN, varCol);
% % plotNURBS(nurbs_temp,[100 100], 1, 1.5*[44 77 32]/255, 1);
% hold on
% % plotNURBS(nurbs_inter,[100 100], 1, 1.5*[44 77 32]/255, 1);
% axis equal
% view(243,14)
% % view(180,90)
% camlight
% axis off
% colorbar
% plotControlPts(nurbs_temp)
% 
% w1 = nurbs_temp.coeffs(1:3,1,3)-nurbs_temp.coeffs(1:3,1,1);
% w2 = nurbs_temp.coeffs(1:3,end,3)-nurbs_temp.coeffs(1:3,end,1);
% nurbsCone = getConeData(abs(-b+h+x2), abs(-b+h+x2) + (L+g2+g3+a-l_ls-58.9)*tan(alpha), L+g2+g3+a-l_ls-58.9, 0, acos(dot(w1,w2)/norm(w1)/norm(w2)), -(L+g2+g3));
% 
% % nurbsCone = elevateNURBSdegree(nurbsCone,[0,1]);
% % nurbsCone = insertKnotsInNURBS(nurbsCone,{[5/12, 13/15] []});
% varCol.colorFun = @(v,xi,eta) log10(uTest_BC(v,xi,eta,nurbs_temp.knots{1},nurbs_temp.knots{2},a,b,L,g2,g3,alpha,[1,numel(unique(Eta)),numel(unique(Eta))]));
% 
% coeffs = nurbsCone.coeffs;
% for i = 1:size(coeffs,2)
%     v1 = coeffs(1:3,i,1);
%     v2 = coeffs(1:3,i,2);
%     
%     g4 = v1(1)+4;
%     
%     coeffs(1:3,i,end) = v1 + (g4-v1(1))/(v2(1)-v1(1))*(v2-v1);
% %     coeffs(1:3,i,2) = (coeffs(1:3,i,1)+coeffs(1:3,i,3))/2;
% end
% nurbsCone.coeffs = coeffs;
% plotNURBS(nurbsCone,[100 100], 1, 1.5*[44 77 32]/255, 1, NaN, varCol);
% w1 = nurbs_temp.coeffs(1:3,1,3)-nurbs_temp.coeffs(1:3,1,1);
% w2 = nurbs_temp.coeffs(1:3,end,3)-nurbs_temp.coeffs(1:3,end,1);
% nurbsCone2 = getConeData(abs(-b+h+x2), abs(-b+h+x2) + g4*tan(alpha), g4, 0, acos(dot(w1,w2)/norm(w1)/norm(w2)), -(L+g2+g3));
% 
% 4
% plotControlPts(nurbsCone,'blue','blue','blue')
%% back rudder
coeffs = ones(4,4,10);
Xi = [0,0,1,2,3,3]/3;
Eta = (nurbs.knots{2}(6:16)-eta_s)/(eta_e-eta_s);
Eta = [0,Eta,1];
coeffs(1:3,1,:) = nurbs.coeffs(1:3,5,5:14);
coeffs(1,2,:) = -delta_s-58.9+a-l_us+l_us*cntrlPts(:,1);
coeffs(3,2,:) = -l_us*cntrlPts(:,2);
coeffs(1,3,:) = -delta_s-58.9+a-l_us+l_us*cntrlPts(:,1);
coeffs(3,3,:) = l_us*cntrlPts(:,2);
coeffs(2,2:3,:) = h_s;
coeffs(1:3,4,:) = nurbs.coeffs(1:3,6,5:14);
nurbs_rudderb = createNURBSobject(coeffs,{Xi, Eta});
nurbs_rudderb = elevateNURBSdegree(nurbs_rudderb,[1,0]);
temp = nurbs_rudderb.coeffs(1,4,end-4:end-1);
nurbs_rudderb.coeffs(1,4,end-4:end-1) = temp(:) - l_us*[0.005; 0.01; 0.023; 0.031];

nurbs_rudderb2 = nurbs_rudderb;
nurbs_rudderb2.coeffs(2,:,:) = -nurbs_rudderb2.coeffs(2,:,:);

nurbs_rudderb3 = nurbs_rudderb;
coeffs = nurbs_rudderb3.coeffs;
for i = 1:size(coeffs,3)
    coeffs(1:3,:,i) = R_x*coeffs(1:3,:,i);
end
nurbs_rudderb3.coeffs = coeffs;
nurbs_rudderb4 = nurbs_rudderb3;
nurbs_rudderb4.coeffs(3,:,:) = -nurbs_rudderb4.coeffs(3,:,:);

%% create remaining hulls for back rudders
coeffs = nurbs.coeffs;
for i = 5:14
    R_x = rotationMatrix(pi/2, 'Xaxis');
    coeffs(1:3,end,i) = R_x*coeffs(1:3,5,i);
    coeffs(1:3,end-1,i) = R_x*coeffs(1:3,4,i);
    R_x = rotationMatrix(-pi/2, 'Xaxis');
    coeffs(1:3,1,i) = R_x*coeffs(1:3,6,i);
    coeffs(1:3,2,i) = R_x*coeffs(1:3,7,i);
end
nurbs.coeffs = coeffs;
% nurbs.coeffs(4,3,5:14) = -nurbs.coeffs(3,6,5:14);

%% Extract upper
Xi = (nurbs.knots{1}(11:end)-xi_r)/(1-xi_r);
Eta = nurbs.knots{2};
nurbs_upper = createNURBSobject(nurbs.coeffs(:,end-3:end,:),{Xi, Eta});

%% Extract inter
Xi = (nurbs.knots{1}(6:13)-0.5)/(xi_r-0.5);
Eta = nurbs.knots{2};
nurbs_inter = createNURBSobject(nurbs.coeffs(:,6:10,:),{Xi, Eta});

%% Extract lower
Xi = nurbs.knots{1}(1:8)/0.5;
Eta = nurbs.knots{2};
nurbs_lower = createNURBSobject(nurbs.coeffs(:,1:5,:),{Xi, Eta});

%% mirror across the xz-plane
nurbs_upper2 = nurbs_upper;
nurbs_upper2.coeffs(2,:,:) = -nurbs_upper2.coeffs(2,:,:);
nurbs_inter2 = nurbs_inter;
nurbs_inter2.coeffs(2,:,:) = -nurbs_inter2.coeffs(2,:,:);
nurbs_lower2 = nurbs_lower;
nurbs_lower2.coeffs(2,:,:) = -nurbs_lower2.coeffs(2,:,:);

nurbsCol = {nurbs_sail,nurbs_rudderb3,nurbs_upper,nurbs_rudder,nurbs_inter,nurbs_rudderb,nurbs_lower,nurbs_rudderb4,nurbs_lower2,nurbs_rudderb2,nurbs_inter2,nurbs_rudder2,nurbs_upper2};


hyp = NaN;
