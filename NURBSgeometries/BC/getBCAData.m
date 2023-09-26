function nurbsColFinal = getBCAData(p,useRefinements)
if nargin < 2
    useRefinements = true;
end
alpha = NaN; % avoid bug
beta = NaN; % avoid bug
container = setBCParameters;
fieldNames = fieldnames(container);
for i = 1:numel(fieldNames)
    eval([fieldNames{i} ' = container.' fieldNames{i} ';']);
end
h = g2*tan(alpha/2); % = g2/tan((pi-alpha)/2)
x2 = g3*tan(alpha);

dx = 58.9-a-L-g2-g3/2;
dy = 0.2;
dy2 = 0.2;

[g,dg,~,dSdxi,dSdeta] = getMainRudderFunctions(b,g2,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m);
xi1 = 0.184706062558860;
s2_r = 1/8;
[~,~,arc] = NACAarcLength2(1,g,dg,dSdxi,dSdeta);
xi2 = invertNACA2_2(s2_r,l_lm,delta_m,x_m,g,dg);
totLength = arc(xi2);
s1 = NACAarcLength2(xi1,g,dg,dSdxi,dSdeta)*s2_r/totLength;
% Xisail = [zeros(1,p+1),xi_t*ones(1,p),(xi_t+st)/2,st*ones(1,p),(st+0.3)/2,0.3*ones(1,p),linspace2(0.3,1,8),ones(1,p+1)];

rudders = cell(1,4);
for rudder_i = 1:2
    if rudder_i == 1 % The top rudder is special, as some extra knots are repeated
        Xi = [zeros(1,p+1),copyVector(linspace2(0,s1,2),p,2).',s1*ones(1,p),(s1+s2_r)/2,s2_r*ones(1,p),linspace2(s2_r,1,4),ones(1,p+1)];
        idx = 4;
    else
        Xi = [zeros(1,p+1),linspace2(0,s1,2),s1*ones(1,p),(s1+s2_r)/2,s2_r*ones(1,p),linspace2(s2_r,1,4),ones(1,p+1)];
        idx = 2;
    end
    [~, I] = getRepeatedKnots(Xi,p);
    
    % Fix front (arc) part of cone part 
    z = g3/2*tan(alpha) + abs(-b+h+x2);
    x = -L-g2-g3/2;
    if rudder_i == 1
        XiTemp = [0,0,0,1,1,2,2,3,3,3]/3;
    else
        XiTemp = [0,0,0,1,1,1];
    end
    n = numel(XiTemp)-(2+1);
    controlPtsTemp = parmArc(XiTemp,30*pi/180);
    controlPts = [x*ones(1,n); z*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
    R_x = rotationMatrix(60*pi/180, 'Xaxis');
    controlPts(1:3,:) = R_x*controlPts(1:3,:);
    nurbsTemp = createNURBSobject(controlPts, XiTemp);
    nurbsTemp = elevateNURBSdegree(nurbsTemp,p-2);
    if rudder_i ~= 1
        nurbsTemp = insertKnotsInNURBS(nurbsTemp,{[1,2]/3});
    end
    
    r_o = @(xi) getBeTSSiConePart(xi,ones(size(xi)),b,g2,g3,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m,xi1,s1,s2_r,'rudder',rudder_i);
%     r_u = @(xi) getBeTSSiConePart(xi,zeros(size(xi)),b,g2,g3,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m,xi1,s1,s2_r,'rudder',rudder_i);
    r_u = @(xi) getBeTSSiConePart(xi,zeros(size(xi)),b,g2,g3,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m,xi1,s1,s2_r,'curve',rudder_i);
    nurbs = addBeTSSiFoil3(l_lm, b_lm, l_um, b_um, x_m, 0, h_m, dx, [dy, dy2], Xi, idx, NaN, s2_r, r_o, r_u);
    topRudder = nurbs;
    nurbs{1}.coeffs = nurbs{1}.coeffs(:,:,1:2);
    nurbs{1}.knots{2} = [0,0,1,1];
    nurbs{1}.number(2) = 2;
    nurbs = elevateNURBSdegree(nurbs,[0,p-1]);
    noNewEtaKnots = 2;
    nurbs = insertKnotsInNURBS(nurbs,{{[] insertUniform2(nurbs{1}.knots{2}, noNewEtaKnots)}});

    nurbs{1}.coeffs(:,1:nurbsTemp{1}.number,1) = nurbsTemp{1}.coeffs(:,end:-1:1);
    
    % Fix back (arc) part of cone part
    x = a-58.9-l_lm;
    dx = L+g2+g3+x;
    z = abs(-b+h+x2)+dx*tan(alpha);
    XiTemp = [0,0,0,1,1,1];
    controlPtsTemp = parmArc(XiTemp,30*pi/180);
    controlPts = [x*ones(1,3); z*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
    R_x = rotationMatrix(60*pi/180, 'Xaxis');
    controlPts(1:3,:) = R_x*controlPts(1:3,:);
    nurbsTemp = createNURBSobject(controlPts, XiTemp);
    nurbsTemp = elevateNURBSdegree(nurbsTemp,p-2);
    nurbsTemp = insertKnotsInNURBS(nurbsTemp,{{insertUniform2(nurbsTemp{1}.knots{1}, noNewEtaKnots)}});
    nurbs{1}.coeffs(:,end,:) = nurbsTemp{1}.coeffs;

    S = @(xi,eta) getBeTSSiConePart(xi,eta,b,g2,g3,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m,xi1,s1,s2_r,'cone',rudder_i);
    P1 = S(s1,1);
    P2 = S(s2_r,1);
    P3 = S(1,1);
    Xitemp = [zeros(1,p+1), (Xi(and(Xi > s1+10*eps,Xi < s2_r-10*eps))-s1)/(s2_r-s1), ones(1,p+1)];
    xiT = aveknt(Xitemp, p+1).';
    nurbs{1}.coeffs(1:3,I(idx):I(idx+1),1) = (repmat(P1,numel(xiT),1) + xiT*(P2-P1)).';
    Xitemp = [zeros(1,p+1), (Xi(Xi > s2_r+10*eps)-s2_r)/(1-s2_r)];
    xiT = aveknt(Xitemp, p+1).';
    nurbs{1}.coeffs(1:3,I(idx+1):end,1) = (repmat(P2,numel(xiT),1) + xiT*(P3-P2)).';

    m = nurbs{1}.number(2);
    for i = 2:I(idx)-1 % fix weigh_ss in eta dir for front part
        weigh_ssTemp = linspace(nurbs{1}.coeffs(4,i,1),1,m);
        nurbs{1}.coeffs(4,i,2:m-1) =  weigh_ssTemp(2:m-1);
    end
    for j = 2:m-1 % fix weigh_ss in xi dir for rear part
        weigh_ssTemp = fliplr(linspace(nurbs{1}.coeffs(4,end,j),1,I(end)-I(idx+1)+1));
        nurbs{1}.coeffs(4,I(idx+1):end-1,j) =  nurbs{1}.coeffs(4,I(idx+1):end-1,j).*weigh_ssTemp(1:end-1);
    end
    
    P1 = S(0,1);
    P2 = S(0,0);
    m = nurbs{1}.number(2);
    etaT = aveknt(nurbs{1}.knots{2}, nurbs{1}.degree(2)+1).';
    nurbs{1}.coeffs(1:3,1,:) = (repmat(P1,m,1)+etaT*(P2-P1)).';

    S = @(xi,eta) getBeTSSiConePart(xi,1-eta,b,g2,g3,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m,xi1,s1,s2_r,'cone',rudder_i);
    if 1
        xiT = aveknt(nurbs{1}.knots{1}, nurbs{1}.degree(1)+1).';
        etaT = aveknt(nurbs{1}.knots{2}, nurbs{1}.degree(2)+1).';
        nurbs = interpolateSurface(nurbs,xiT,etaT,S);
    else
        nurbs = leastSquares2D(nurbs,S,3);
    end
    nurbsMir = nurbs;
    nurbsMir{1}.coeffs(2,:,:) = -nurbsMir{1}.coeffs(2,:,:);
    nurbsMir = flipNURBSparametrization(nurbsMir,2);

    topRudder{1}.coeffs = topRudder{1}.coeffs(:,:,2:end-1);
    topRudder{1}.knots{2} = [0,0,1,2,3,4,4]/4;
    topRudder{1}.number(2) = 5;
    topRudder{1}.coeffs(:,:,end) = nurbsMir{1}.coeffs(:,:,1);
    topRudder{1}.coeffs(:,:,1) = nurbs{1}.coeffs(:,:,end);


    topRudder = elevateNURBSdegree(topRudder,[0,p-1]);
    if useRefinements
	    topRudder = insertKnotsInNURBS(topRudder,{{[] [linspace2(0,1,4), linspace2(2,3,4)]/3}});
    end
    
    rudders(rudder_i) = glueNURBS([nurbs,topRudder,nurbsMir],2);
end
rudders(2) = rotateNURBS(rudders(2),'theta',pi/2,'rotAxis', [1,0,0]);
rudders(3) = rotateNURBS(rudders(2),'theta',pi/2,'rotAxis', [1,0,0]);
rudders(4) = rotateNURBS(rudders(3),'theta',pi/2,'rotAxis', [1,0,0]);

b2 = b-h;
x_b3 = -L-g2-g3/2;
b3 = (x_b3+L+g2+cot(alpha)*(b-h))/cot(alpha);
x_b4 = a-58.9-l_lm;
b4 = (x_b4+L+g2+cot(alpha)*(b-h))/cot(alpha);
b5 = abs(-b+h+x2);

XiTemp = [0,0,0,1,1,2,2,3,3,3]/3;
n = numel(XiTemp)-(2+1);
coeffs = zeros(4,n,2);
controlPtsTemp = parmArc(XiTemp,30*pi/180);
coeffs(:,:,1) = [x_b4*ones(1,n); b4*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
coeffs(:,:,2) = [x_b3*ones(1,n); b3*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];

conePart = cell(1,4);
conePart(1) = createNURBSobject(coeffs,{XiTemp, [0,0,1,1]});
conePart(1) = insertKnotsInNURBS(conePart(1),{[] 1-(s2_r-s1)/(1-s1)});
conePart(1) = elevateNURBSdegree(conePart(1),[p-2,p-1]);

XiTemp = [0,0,0,1,1,1];
n = numel(XiTemp)-(2+1);
coeffs = zeros(4,n,2);
controlPtsTemp = parmArc(XiTemp,30*pi/180);
coeffs(:,:,1) = [x_b4*ones(1,n); b4*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
coeffs(:,:,2) = [x_b3*ones(1,n); b3*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];

conePart(2) = createNURBSobject(coeffs,{XiTemp, [0,0,1,1]});
conePart(2) = insertKnotsInNURBS(conePart(2),{[] 1-(s2_r-s1)/(1-s1)});
conePart(2) = elevateNURBSdegree(conePart(2),[p-2,p-1]);
if useRefinements
    conePart(2) = insertKnotsInNURBS(conePart(2),{[1,2]/3 []});
end


conePart(1) = rotateNURBS(conePart(1),'theta',pi/2+30*pi/180,'rotAxis', [1,0,0]);
conePart(2) = rotateNURBS(conePart(2),'theta',pi+30*pi/180,'rotAxis', [1,0,0]);
conePart(3) = rotateNURBS(conePart(2),'theta',pi/2,'rotAxis', [1,0,0]);
conePart(4) = rotateNURBS(conePart(1),'theta',-pi/2,'rotAxis', [1,0,0]);
XiTemp = rudders{3}.knots{1};
XiTemp = [zeros(1,p+1), 1-fliplr((XiTemp(and(XiTemp > s1+eps,XiTemp < 1))-s1)/(1-s1)),ones(1,p+1)];
for i = 1:4
    conePart(i) = insertKnotsInNURBS(conePart(i),{[] [linspace2(0,1-(s2_r-s1)/(1-s1),4), 1-(s1+s2_r)/2]});
    conePart{i}.knots{2} = XiTemp;
    for j = 1:size(conePart{i}.coeffs,2)
        x = rudders{3}.coeffs(1,end:-1:I(idx),end);
        x = reshape(x,1,1,numel(x));
        x_old = conePart{i}.coeffs(1,j,:);
        conePart{i}.coeffs(1,j,:) = x;
        b_x = (x+L+g2+cot(alpha)*(b-h))/cot(alpha);
        b_xold = (x_old+L+g2+cot(alpha)*(b-h))/cot(alpha);
        conePart{i}.coeffs(2:3,j,:) = conePart{i}.coeffs(2:3,j,:).*repmat(b_x./b_xold,2,1);
    end
end
    

XiTemp = [0,0,0,3,3,4,4,5,5,6,6,9,9,12,12,13,13,14,14,15,15,18,18,21,21,24,24,27,27,30,30,33,33,36,36,36]/36;
n = numel(XiTemp)-(2+1);
coeffs = zeros(4,n,2);
controlPtsTemp = parmArc(XiTemp,2*pi);
x = -L-g2-g3;
controlPts = [x*ones(1,n); b5*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
coeffs(:,:,2) = controlPts;
coeffs(1,:,1) = x;
coeffs(4,:,1) = coeffs(4,:,2);
backDisk = createNURBSobject(coeffs,{XiTemp, [0,0,1,1]});
backDisk = elevateNURBSdegree(backDisk,[p-2,p-1]);
if useRefinements
    refineXi = [1,2,7,8,10,11,16,17,19,20,22,23,25,26,28,29,31,32,34,35]/36;
else
    refineXi = [1,2,7,8,10,11,16,17,19,20,25,26,28,29,34,35]/36;
end
backDisk = insertKnotsInNURBS(backDisk,{refineXi []});

x = a-58.9-l_lm;
coeffs(:,:,1) = coeffs(:,:,2);
controlPts = [x*ones(1,n); b4*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
coeffs(:,:,2) = controlPts;
backCone = createNURBSobject(coeffs,{XiTemp, [0,0,1,1]});
backCone = elevateNURBSdegree(backCone,[p-2,p-1]);
backCone = insertKnotsInNURBS(backCone,{refineXi []});

XiTemp = [0,0,0,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,18,18,21,21,24,24,27,27,30,30,33,33,36,36,36]/36;
n = numel(XiTemp)-(2+1);
coeffs = zeros(4,n,2);
controlPtsTemp = parmArc(XiTemp,2*pi);
x_b3 = -L-g2-g3/2;
coeffs(:,:,1) = [x_b3*ones(1,n); b3*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
x_b2 = -L-g2;
coeffs(:,:,2) = [x_b2*ones(1,n); b2*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
frontCone = createNURBSobject(coeffs,{XiTemp, [0,0,1,1]});
frontCone = elevateNURBSdegree(frontCone,[p-2,p-1]);
if useRefinements
    refineXi = [1,2,16,17,19,20,22,23,25,26,28,29,31,32,34,35]/36;
    refineEta = linspace2(0,1,4);
else
    refineXi = [1,2,16,17,19,20,25,26,28,29,34,35]/36;
    refineEta = [];
end
frontCone = insertKnotsInNURBS(frontCone,{refineXi refineEta});

nurbsColFinal = cell(1,36);
nurbsColFinal(1) = backDisk;
nurbsColFinal(2) = backCone;
nurbsColFinal(3:2:9) = rudders;
nurbsColFinal(4:2:10) = conePart;
nurbsColFinal(11) = frontCone;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nurbsCol = cell(1,23);
dx_rudder = 0.05;
repLoc = zeros(1,8);
s2_dr = 1/8;
dx_sail = 0.2;
s2_s = 1/8;
repLoc(1:4) = [dx_rudder+x_d, x_d-s2_dr*l_ld, x_d-0.3*l_ld,x_d-l_ld];
repLoc(5:8) = [dx_sail+x_s, x_s-s2_s*l_ls, x_s-0.3*l_ls, x_s-l_ls];

repLoc = (L+fliplr(repLoc))/L;
nurbsCol(1) = getBC_modData2(a, b, L, g2, g3, alpha, beta, c, s,repLoc,true);
nurbsCol(1) = elevateNURBSdegree(nurbsCol(1),[p-2,p-2]);
nurbsCol2 = explodeNURBS(nurbsCol(1),2);
nurbsCol3 = glueNURBS(nurbsCol2(2:10),2);
nurbsCol3 = explodeNURBS(nurbsCol3,1);
nurbsCol(1:3) = nurbsCol3(1:3);
nurbsCol(13:end) = nurbsCol3(10:end);
nurbsColFinal(12) = nurbsCol2(1);
nurbsColFinal(end) = nurbsCol2(end);

%% add sail
xi1 = 0.155350810458409;
xi2 = sqrt(s2_s);
t = b_ls/l_ls;
totLength = NACAarcLength(xi2,t);
s1 = NACAarcLength(xi1,t)*s2_s/totLength;
idx = 2;

deck = glueNURBS(nurbsCol3(6:7),1);
deck = explodeNURBS(deck,2);
nurbsCol(7) = deck(1);
nurbsCol(9) = glueNURBS(deck(5:9),2);
Xisail = [zeros(1,p+1),s1*ones(1,p),(s1+s2_s)/2,s2_s*ones(1,p),(s2_s+0.3)/2,0.3*ones(1,p),linspace2(0.3,1,8),ones(1,p+1)];

r_o = @(xi) getBeTSSiSail(xi,ones(size(xi)),a,c,l_ls,l_us,b_ls,b_us,h_s,delta_s, s2_s);
r_u = @(xi) getBeTSSiSail(xi,zeros(size(xi)),a,c,l_ls,l_us,b_ls,b_us,h_s,delta_s, s2_s);
sail = addBeTSSiFoil3(l_ls, b_ls, l_us, b_us, x_s, c, h_s, dx_sail, [0,0], Xisail, idx, NaN, s2_s, r_o, r_u);
sail = elevateNURBSdegree(sail,[0,p-1]);
nurbsCol(8) = sail;

%% add depthrudders
y0 = s+dyPanels;
y1 = s+2*dyPanels;

z0 = P_panels(y0);
z1 = P_panels(y1);
idx = 2;

xi1 = 0.201905012203230;
[~,~,arc] = NACAarcLength2_depth(1,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,'rudderupper');
xi2 = invertNACA2_depth_2(s2_dr,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,'rudderupper');
% xi2 =   invertNACA2_depth(s2_dr,a,b,l_ls,l_us,b_ls,b_us,h_s,delta_o,beta,s,c,'rudderupper',s2_dr);
totLength = arc(xi2);
s1 = NACAarcLength2_depth(xi1,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,'rudderupper')*s2_dr/totLength;

Xirudder = [zeros(1,p+1),s1*ones(1,p),s2_dr*ones(1,p),0.3*ones(1,p),linspace2(0.3,1,3),ones(1,p+1)];
r_o1 = @(xi) getBeTSSiPanelPart(xi,ones(size(xi)),b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2_dr,'rudderupper','curveupper');
r_o2 = @(xi) getBeTSSiPanelPart(xi,ones(size(xi)),b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2_dr,'rudderupper','curveupper');
r_u = @(xi) getBeTSSiPanelPart(xi,zeros(size(xi)),b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2_dr,'rudderupper','curveupper');
leftDepthRudder = addBeTSSiFoil3(l_ld, b_ld, l_ud, b_ud, x_d, s, h_d, dx_rudder,[0,c-z1-b_ld], Xirudder, idx,NaN,s2_dr, {r_o1,r_o2}, r_u, 3, c-b_ld/2);

r = @(xi) getBeTSSiPanelPart(xi,zeros(size(xi)),b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2_dr,'curveupper');
[xi_rep, I] = getRepeatedKnots(Xirudder,p);

XiTemp = [zeros(1,p+1), ones(1,p+1)];
xiT1 = aveknt(XiTemp, p+1);
X1 = @(xi) x_d - s2_dr*l_ld - xi*(0.3-s2_dr)*l_ld;
XiTemp2 = [zeros(1,p+1), (Xirudder(Xirudder > 0.3+10*eps)-0.3)/(1-0.3)];
xiT2 = aveknt(XiTemp2, p+1);
X2 = @(xi) x_d - 0.3*l_ld - xi*(1-0.3)*l_ld;
temp1 = X1(xiT1);
temp2 = X2(xiT2);

dofsToRemove{1} = sort([I,I(idx+2)-1,I(idx+2)+1]);
dofsToRemove{2} = sort([I,I(idx+2)-1,I(idx+2)+1]);
dofsToRemove{3} = sort([I,I(idx+2)-1,I(idx+2)+1]);
temp = r(xi_rep.');
values{1} = [temp(1:idx+1,1); temp1(end-1); temp(idx+2,1); temp2(2); temp(idx+3:end,1)];
% values{1} = [temp(1,1); temp(1:idx+1,1); temp1(end-1); temp(idx+2,1); temp2(2); temp(idx+3:end,1)];
values{2} = temp(:,2);
values{2} = [values{2}(1:idx+1); s*ones(3,1); values{2}(idx+3:end)];
values{3} = temp(:,3);
values{3} = [values{3}(1:idx+1); c*ones(3,1); values{3}(idx+3:end)];
nurbs_upper = leastSquares1D(Xirudder,p,r,3,ones(1,I(end)),dofsToRemove,values);
r = @(xi) getBeTSSiPanelPart(xi,zeros(size(xi)),b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2_dr,'curvelower','curveupper');
nurbs_lower = leastSquares1D(Xirudder,p,r,3);

coeffs = leftDepthRudder{1}.coeffs;
coeffs(2,:,end) = y0 + (y1-y0)/(z1-z0)*(coeffs(3,:,end)-z0);
coeffs(1:3,:,end-1) = nurbs_lower{1}.coeffs(1:3,:);
y0 = s;
y1 = s+dyPanels;
z0 = P_panels(y0);
z1 = P_panels(y1);
coeffs(2,:,1) = y0 + (y1-y0)/(z1-z0)*(coeffs(3,:,1)-z0);
coeffs(1:3,:,2) = nurbs_upper{1}.coeffs(1:3,:);
leftDepthRudder{1}.coeffs = coeffs;

sidePanel = glueNURBS(nurbsCol3(4:5),1);
sidePanel = explodeNURBS(sidePanel,2);
nurbsCol(4) = glueNURBS(sidePanel(1:5),2);
nurbsCol(6) = sidePanel(9);
leftDepthRudder = elevateNURBSdegree(leftDepthRudder,[0,p-1]);
leftDepthRudder = flipNURBSparametrization(leftDepthRudder,2);
nurbsCol(5) = leftDepthRudder;

nurbsCol(10) = flipNURBSparametrization(mirrorNURBS(nurbsCol(4),'y'),1);
nurbsCol(11) = flipNURBSparametrization(mirrorNURBS(nurbsCol(5),'y'),2);
nurbsCol(12) = flipNURBSparametrization(mirrorNURBS(nurbsCol(6),'y'),1);

nurbsColFinal(13:end-1) = nurbsCol;

nurbsColFinal = makeUniformNURBSDegree(nurbsColFinal,p);

%% Refine to reduce aspect ratios and to make patches compatible

% refine side panels and lower part
if useRefinements
    refineXi = [linspace2(0,1/9,14),  linspace2(1/9,2/9,8), linspace2(2/9,3/9,1), ...
                linspace2(3/9,4/9,1), linspace2(4/9,5/9,4), linspace2(5/9,6/9,3), linspace2(8/9,1,3)];
else
    refineXi = [linspace2(1/9,2/9,8), linspace2(2/9,3/9,1),linspace2(3/9,4/9,1), linspace2(5/9,6/9,3)];
end
for i = [13:15, 25:35]
    nurbsColFinal(i) = insertKnotsInNURBS(nurbsColFinal(i),{[] refineXi});
end

% refine upper panels
if useRefinements
    refineXi = [linspace2(0,1/5,14),linspace2(1/5,2/5,8), linspace2(2/5,3/5,1), linspace2(3/5,4/5,1), linspace2(4/5,1,4)];
else
    refineXi = [linspace2(1/5,2/5,8), linspace2(2/5,3/5,1), linspace2(3/5,4/5,1)];
end
for i = [16, 22]
    nurbsColFinal(i) = insertKnotsInNURBS(nurbsColFinal(i),{[] refineXi});
end

if useRefinements
    % refine depth rudders
    for i = [17,23]
        nurbsColFinal(i) = insertKnotsInNURBS(nurbsColFinal(i),{[] [linspace2(1/5,2/5,3), linspace2(3/5,4/5,3)]});
    end
    % refine rear deck
    nurbsColFinal(19) = insertKnotsInNURBS(nurbsColFinal(19),{[] linspace2(0,1,14)});
    % refine upper front side panels
    for i = [18, 24]
        nurbsColFinal(i) = insertKnotsInNURBS(nurbsColFinal(i),{[] linspace2(0,1,3)});
    end
    % refine sail
    nurbsColFinal(20) = insertKnotsInNURBS(nurbsColFinal(20),{[] [linspace2(1/5,2/5,3), linspace2(3/5,4/5,3)]});
end

% refine front deck
if useRefinements
    refineXi = [linspace2(0,1/5,4), linspace2(1/5,2/5,3), linspace2(4/5,1,3)];
else
    refineXi = linspace2(1/5,2/5,3);
end
nurbsColFinal(21) = insertKnotsInNURBS(nurbsColFinal(21),{[] refineXi});

% refine in xi-direction for lower part
if useRefinements
    indices = 28:35;
else
    indices = [28,29,31,32,34,35];
end
for i = indices
    nurbsColFinal(i) = insertKnotsInNURBS(nurbsColFinal(i),{linspace2(0,1,2) []});
end

% refine transition part
if useRefinements
    Xi_temp = nurbsColFinal{12}.knots{1};
    refineXi = insertUniform2([1/3, Xi_temp(Xi_temp>=1/3)],2);
    refineEta = [linspace2(0,0.5,3), linspace2(0.5,1,3)];
else
    refineXi = [13,14,16,17,22,23,25,26,31,32,34,35]/36;
    refineEta = [];
end
nurbsColFinal(12) = insertKnotsInNURBS(nurbsColFinal(12),{refineXi refineEta});

% refine bow
if useRefinements
    Xi_temp = nurbsColFinal{end}.knots{1};
    refineXi = insertUniform2([1/3, Xi_temp(Xi_temp>=1/3)],2);
    refineEta = linspace2(0,1,14);
else
    refineXi = [13,14,16,17,22,23,25,26,31,32,34,35]/36;
    refineEta = [];
end

nurbsColFinal(end) = insertKnotsInNURBS(nurbsColFinal(end),{refineXi refineEta});

nurbsColFinal = explodeNURBS(nurbsColFinal,2);
nurbsColFinal = explodeNURBS(nurbsColFinal,1);
nurbsColFinal = explodeNURBS(nurbsColFinal,2);
nurbsColFinal = explodeNURBS(nurbsColFinal,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nurbs = interpolateSurface(nurbs,xiT,etaT,f)

P = nurbs{1}.coeffs;
n_xi = nurbs{1}.number(1);
n_eta = nurbs{1}.number(2);

p_xi = nurbs{1}.degree(1);
p_eta = nurbs{1}.degree(2);
degree = nurbs{1}.degree;

t_xi = nurbs{1}.knots{1};
t_eta = nurbs{1}.knots{2};
knots = nurbs{1}.knots;

uniqueXi = unique(t_xi);
for i = 1:n_xi
    I = find(abs(uniqueXi - xiT(i)) < 10*eps);
    if ~isempty(I)
        xiT(i) = uniqueXi(I);
    end
end
uniqueEta = unique(t_eta);
for j = 1:n_eta
    I = find(abs(uniqueEta - etaT(j)) < 10*eps);
    if ~isempty(I)
        etaT(j) = uniqueEta(I);
    end
end

% Create interpolation points
Q = P(1:3,:,:);
xi_arr = aveknt(t_xi, p_xi+1);
eta_arr = aveknt(t_eta, p_eta+1);
[XI,ETA] = ndgrid(xi_arr(2:end-1),eta_arr(2:end-1));
Q(:,2:end-1,2:end-1) = reshape(f(XI(:),ETA(:)).',3,n_xi-2,n_eta-2);

i_xi_arr = zeros(n_xi,1);
for i = 1:n_xi
    i_xi_p = findKnotSpan(n_xi, p_xi, xiT(i), t_xi);
    i_xi_arr(i) = i_xi_p;
end

j_eta_arr = zeros(n_eta,1);
for j = 1:n_eta
    j_eta_p = findKnotSpan(n_eta, p_eta, etaT(j), t_eta);
    j_eta_arr(j) = j_eta_p;
end

index = zeros(n_xi*n_eta,2);
counter = 1;
for j = 1:n_eta
    for i = 1:n_xi   
        index(counter,:) = [i, j];
        counter = counter + 1;
    end
end
n_en = (p_xi+1)*(p_eta+1);
weigh_ss = reshape(nurbs{1}.coeffs(4,:,:),n_xi*n_eta,1);
values = zeros(n_xi*n_eta,n_en);
n_idx = zeros(n_xi*n_eta,n_en);
m_idx = zeros(n_xi*n_eta,n_en);
for a = 1:n_xi*n_eta
    i = index(a,1);
    j = index(a,2);

    i_xi = i_xi_arr(i);
    j_eta = j_eta_arr(j);
    m = i + (j-1)*n_xi;

    temp = zeros(1,n_en);
    counter = 1;
    for jj = 1:p_eta+1
        j_tilde = j_eta - p_eta + jj - 1;
        n = i_xi - p_xi + (1:p_xi+1) - 1 + (j_tilde-1)*n_xi;
        
        temp(counter:counter+p_xi) = n;
        counter = counter + p_xi+1;
    end
    xi = [xiT(i), etaT(j)];
    I = findKnotSpans(degree, xi, knots);
    R = NURBSbasis(I, xi, degree, knots, weigh_ss(temp));
    values(a,:) = R{1};
    n_idx(a,:) = temp;
    m_idx(a,:) = repmat(m,1,n_en);
end
A = sparse(m_idx,n_idx,values);
indices = sort([2:n_xi-1, 1:n_xi:n_xi*n_eta, n_xi:n_xi:n_xi*n_eta, n_xi*(n_eta-1)+2:n_xi*n_eta-1]);
A(indices,:) = [];
P = reshape(P,4,n_xi*n_eta);

Q = reshape(Q,3,n_xi*n_eta);
Q(:,indices) = [];
Q = Q - (A(:,indices)*P(1:3,indices).').';
A(:,indices) = [];
P_tilde = zeros(size(Q));
P_tilde(1,:) = (A\Q(1,:)')';
P_tilde(2,:) = (A\Q(2,:)')';
P_tilde(3,:) = (A\Q(3,:)')';
P(1:3,setdiff(1:n_xi*n_eta, indices')) = P_tilde;
nurbs{1}.coeffs = reshape(P,4,n_xi,n_eta);
