function nurbs = getBC_modData2(a, b, L, g2, g3, alpha, beta, c, s, repLoc)
h = g2*tan(alpha/2); % = g2/tan((pi-alpha)/2)
t1 = 2*pi-beta;
npts = 6;
p2 = t1/(2*npts);

x1 = tan(alpha/2)*(g2/tan(alpha) + h);
lowerCpts = [  -L-g2,       b-h,        1;
               -L-x1,       b,          cos(alpha/2);
               -L,          b,          1;
               -L/2,        b,          1;
                0,          b,          1;
                a,          b,          1/sqrt(2);
                a,          0,          1].';

Xi1 = [0,0,0,3,3,6,6,9,9,12,12,15,15,18,18,21,21,24,24,24]/24;
Eta = [0 0 0 1 1 2 2 3 3 3]/3;
ctrlPtsXi = parmArc(Xi1,beta);

coeffs = calcTensorRotCtrlPts(ctrlPtsXi,lowerCpts);
nurbs_lower = rotateNURBS(createNURBSobject(coeffs,{Xi1, Eta}),pi-3*p2,'Xaxis');
nurbs_lower = insertKnotsInNURBS(nurbs_lower{1},{[] 1/6});

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
backTop = rotateNURBS(createNURBSobject(coeffs,{Xi2, [0,0,1,1]}),30*pi/180,'Xaxis');
backTop = backTop{1};

%% Upper main body
C_4 = c + b*cos(beta/2);
C_3 = b*sin(beta/2)-s;
C_2 = (2*C_4+C_3*tan(beta/2))/C_3^3;
C_1 = -(3*C_4+C_3*tan(beta/2))/C_3^2;

P = @(y) c + C_1*(y-s).^2 + C_2*(y-s).^3;

dy = C_3/(npts-1);

coeffs_mainTop = zeros(4,4*npts+1,5);

y = s + (npts-1)*dy;
for i = 1:npts
    line = zeros(5,4);
    line(1,:)   = [-L,   y, P(y), 1];
    line(2,:)   = [-L/2, y, P(y), 1];
    line(3:5,:) = [0,    y, P(y), 1;
                   a,    y, P(y), 1/sqrt(2);
                   a,    0,   0 , 1];
    coeffs_mainTop(:,2*i-1,:) = line.';
    y = y - dy;
end
coeffs_mainTop(:,2*npts+1,:) = coeffs_mainTop(:,2*npts-1,:);
coeffs_mainTop(2,2*npts+1,:) = 0;

for i = 1:npts
    coeffs_mainTop(1:3,2*i,:) = (coeffs_mainTop(1:3,2*i-1,:)+coeffs_mainTop(1:3,2*i+1,:))/2;
    coeffs_mainTop(4,2*i,:) = coeffs_mainTop(4,2*i-1,:);
end

coeffs_mainTop(:,2*npts+2:end,:) = coeffs_mainTop(:,2*npts:-1:1,:);
coeffs_mainTop(2,2*npts+2:end,:) = -coeffs_mainTop(2,2*npts+2:end,:);

%% Create transition and combine into upper main body
coeffs_top = zeros(4,25,8);
coeffs_top(:,:,4:end) = coeffs_mainTop;
coeffs_top(:,:,1) = backTop.coeffs(:,:,2);
w2 = cos(p2/2);
weigh_ss = linspace(w2,1,4);
for i = 1:2*npts+1
    pts = zeros(2,4);
    pts(2,1:3) = coeffs_top(1:3,i,4);
    if mod(i,2)
        pts(1,4) = 1;
        pts(2,4) = 1;
    else
        pts(1,4) = weigh_ss(2);
        pts(2,4) = weigh_ss(3);
    end
    pts(2,1) = nurbs_lower.coeffs(1,1,3);
    
    v1 = backTop.coeffs(1:3,i,1);
    v2 = backTop.coeffs(1:3,i,2);
    
    g4 = nurbs_lower.coeffs(1,1,2)-nurbs_lower.coeffs(1,1,1);
    
    pts(1,1:3) = (v2-v1)/norm(v2-v1)*g4/cos(alpha)+v2;
    coeffs_top(:,i,2:3) = pts.';
end
coeffs_top(4,:,2:3) = coeffs_top(4,:,2:3)*0.5*(1+cos(alpha/2));
coeffs_top(:,(2*npts+2):end,2:3) = coeffs_top(:,2*npts:-1:1,2:3);
coeffs_top(2,(2*npts+2):end,2:3) = -coeffs_top(2,(2*npts+2):end,2:3);


nurbs_top = createNURBSobject(coeffs_top,{Xi2, nurbs_lower.knots{2}});
nurbs = glueNURBS({nurbs_top,nurbs_lower},'xi');
nurbs.knots{1} = [Xi2(1:end-1)/3,1/3+Xi1(4:end)*2/3]; % make knot vector based on the arclength of the rotationally summetric back of submarine

if nargin == 10
    nurbs = insertKnotsInNURBS(nurbs,{[] copyVector(1/3+repLoc/3,2,1)}); % add degrees of freedom to enable C0 lines at rudders and sail
end
