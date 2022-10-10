function nurbsCol = getBC_modData(a, b, L, g2, g3, alpha, beta, c, s)

h = g2*tan(alpha/2); % = g2/tan((pi-alpha)/2)
x2 = g3*tan(alpha);
t1 = 2*pi-beta;
npts = 6;
p2 = t1/(2*npts);
b2 = b-h;
b5 = abs(-b+h+x2);

%% Create back disk and cone
XiTemp = [0,0,0,1,1,2,2,3,3,3]/3;
n = numel(XiTemp)-(2+1);
coeffs = zeros(4,n,2);
controlPtsTemp = parmArc(XiTemp,2*pi);
x = -L-g2-g3;
controlPts = [x*ones(1,n); b5*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
coeffs(:,:,2) = controlPts;
coeffs(1,:,1) = x;
coeffs(4,:,1) = coeffs(4,:,2);
backDisk = rotateNURBS(createNURBSobject(coeffs,{XiTemp, [0,0,1,1]}), 'theta', 3*p2, 'rotAxis', 'Xaxis');

x = -L-g2;
coeffs(:,:,1) = coeffs(:,:,2);
controlPts = [x*ones(1,n); b2*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
coeffs(:,:,2) = controlPts;
backCone = rotateNURBS(createNURBSobject(coeffs,{XiTemp, [0,0,1,1]}), 'theta', 3*p2, 'rotAxis', 'Xaxis');
nurbsCol(1) = glueNURBS({backDisk{1}, backCone{1}}, 2);

%% Create lower curve
x1 = tan(alpha/2)*(g2/tan(alpha) + h);
lowerCpts = [  -L-g2,       b-h,        1;
               -L-x1,       b,          cos(alpha/2);
               -L,          b,          1;
               -L/2,        b,          1;
                0,          b,          1;
                a,          b,          1/sqrt(2);
                a,          0,          1].';

Xi1 = [0,0,0,1,1,2,2,2]/2;
Eta = [0 0 0 1 1 2 2 3 3 3]/3;
% ctrlPtsXi = parmArc(Xi1,beta);
nurbs = createNURBSobject(lowerCpts,Eta);
nurbs = revolveNURBS(nurbs, 'Xi',Xi1,'theta',beta,'rotAxis', 'Xaxis');
% coeffs = calcTensorRotCtrlPts(ctrlPtsXi,lowerCpts);
nurbs = rotateNURBS(nurbs,'theta', pi-3*p2,'rotAxis', 'Xaxis');
nurbs = permuteNURBS(nurbs,[2,1]);
nurbs = explodeNURBS(nurbs,2);
nurbsCol(2) = nurbs(1);
nurbsCol(3) = createNURBSobject(nurbs{2}.coeffs(:,:,[1,3]),{Xi1, [0,0,1,1]});
nurbsCol(4) = nurbs(3);

%% Create back top
b2 = b-h;
x_b3 = -L-g2-g3/2;
b3 = (x_b3+L+g2+cot(alpha)*(b-h))/cot(alpha);
Xi2 = [0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,12]/12;
n = numel(Xi2)-(2+1);
coeffs = zeros(4,n,2);
controlPtsTemp = parmArc(Xi2,beta/2);
x_b3 = -L-g2-g3/2;
coeffs(:,:,1) = [x_b3*ones(1,n); b3*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
x_b2 = -L-g2;
coeffs(:,:,2) = [x_b2*ones(1,n); b2*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
backTop = rotateNURBS(createNURBSobject(coeffs,{Xi2, [0,0,1,1]}),'theta', 30*pi/180,'rotAxis', 'Xaxis');
backTop = backTop{1};

%% Upper main body
C_4 = c + b*cos(beta/2);
C_3 = b*sin(beta/2)-s;
C_2 = (2*C_4+C_3*tan(beta/2))/C_3^3;
C_1 = -(3*C_4+C_3*tan(beta/2))/C_3^2;

P = @(y) c + C_1*(y-s).^2 + C_2*(y-s).^3;

dy = C_3/(npts-1);

coeffs_mainTop = zeros(4,2*npts,4);

y = s + (npts-1)*dy;
for i = 1:npts
    line = zeros(4,4);
    line(1,:)   = [-L,   y, P(y), 1];
    line(2:4,:) = [0,    y, P(y), 1;
                   a,    y, P(y), 1/sqrt(2);
                   a,    0,   0 , 1];
    coeffs_mainTop(:,i,:) = line.';
    y = y - dy;
end

coeffs_mainTop(:,npts+1:end,:) = coeffs_mainTop(:,npts:-1:1,:);
coeffs_mainTop(2,npts+1:end,:) = -coeffs_mainTop(2,npts+1:end,:);
Xi3 = [0,0,1,2,3,4,5,7,8,9,10,11,12,12]/12;
nurbs_mainTop = createNURBSobject(coeffs_mainTop(:,:,1:2),{Xi3, [0,0,1,1]});
nurbsCol(5) = nurbs_mainTop;
nurbsCol(6) = createNURBSobject(coeffs_mainTop(:,:,2:4),{Xi3, [0,0,0,1,1,1]});

%% Create transition
coeffs_top = zeros(4,4*npts+1,4);
nurbs_lower = insertKnotsInNURBS(nurbsCol(2),{[] 0.5});
coeffs_top(:,:,1) = backTop.coeffs(:,:,2);
nurbs_mainTop = insertKnotsInNURBS(nurbs_mainTop,{0.5 []});
nurbs_mainTop = elevateNURBSdegree(nurbs_mainTop,[1,0]);
coeffs_top(:,:,end) = nurbs_mainTop{1}.coeffs(:,:,1);
w2 = cos(p2/2);
weigh_ss = linspace(w2,1,4);
for i = 1:2*npts+2
    pts = zeros(2,4);
    pts(2,1:3) = coeffs_top(1:3,i,4);
    if mod(i,2)
        pts(1,4) = 1;
        pts(2,4) = 1;
    else
        pts(1,4) = weigh_ss(2);
        pts(2,4) = weigh_ss(3);
    end
    pts(2,1) = nurbs_lower{1}.coeffs(1,1,3);
    
    v1 = backTop.coeffs(1:3,i,1);
    v2 = backTop.coeffs(1:3,i,2);
    
    g4 = nurbs_lower{1}.coeffs(1,1,2)-nurbs_lower{1}.coeffs(1,1,1);
    
    pts(1,1:3) = (v2-v1)/norm(v2-v1)*g4/cos(alpha)+v2;
    coeffs_top(:,i,2:3) = pts.';
end
coeffs_top(4,:,2:3) = coeffs_top(4,:,2:3)*0.5*(1+cos(alpha/2));
coeffs_top(:,(2*npts+2):end,2:3) = coeffs_top(:,2*npts:-1:1,2:3);
coeffs_top(2,(2*npts+2):end,2:3) = -coeffs_top(2,(2*npts+2):end,2:3);

nurbsCol(7) = createNURBSobject(coeffs_top,{Xi2, [0,0,0,0.5,1,1,1]});
