function nurbs = addBeTSSiFoil2(l_ls, b_ls, l_us, b_us, h_s, delta_s, x_0, alpha,dx,dy,dy2)

t = b_ls/l_ls;
f = @(xi) getNACA(xi,t);
% setBCParameters
% f_old = @(x) 5*t*(a0*sqrt(x)-a1*x-a2*x.^2+a3*x.^3-a4*x.^4);
% 100*sqrt(integral(@(x)(f_old(x)-f(x)).^2,0,1)/integral(@(x)f_old(x).^2,0,1))

% return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_C = [-dx,mean([dy,dy2])+b_ls/2]/l_ls;
P = @(xi) [xi^2,getNACA(xi,t,0)];
dP = @(xi) [2*xi,getNACA(xi,t,1)];
d2P = @(xi) [2,getNACA(xi,t,2)];
maxIter = 100;
Eps = 1e-15;
xi_t = 0.168967783470083;
g = @(xi) dot(dP(xi),P(xi)-P_C);
dgdxi = @(xi) dot(d2P(xi),P(xi)-P_C)+dot(dP(xi),dP(xi));
xi_t = newtonsMethod(g,dgdxi,xi_t,maxIter,Eps);
xi_T = [0,xi_t,sqrt(0.3),1];
N = [1, 5, 9, 16];
p = 2;
[cntrlPts, Xi] = getNACAapprox(t,p,N,xi_T,f(xi_T));

coeffs = ones(3,size(cntrlPts,1));
coeffs(1:2,:) = cntrlPts.';
%%%%%%%%%%%%
% Xi = [zeros(1,p+1), x(6+3*(N-1)), x(6+3*(N-1):end) ones(1,p+1)];
% coeffs = reshape(coeffs,3,numel(coeffs)/3);
% nurbsAppr = createNURBSobject(coeffs,Xi);
% P = @(xi) evaluateNURBS(nurbsAppr, xi);
% axis equal
% 100*sqrt(integral(@(x)(f_approx(sqrt(x),P)-f(x)).^2,0,1,'ArrayValued',true)/integral(@(x)f(x).^2,0,1,'ArrayValued',true))

nurbsAppr2 = createNURBSobject(coeffs,Xi);
coeffs2 = nurbsAppr2.coeffs;
coeffs2(:,:,2) = nurbsAppr2.coeffs;
coeffs2(1,1:N(2),2) = -dx/l_ls;
coeffs2(2,N(2):end,2) = (dy+b_ls/2)/l_ls;
% xi_arr = aveknt(Xi, p+1);
% dP = @(xi) evaluateNURBS_deriv2(nurbsAppr, xi, 'xi');
% d2P = @(xi) evaluateNURBS_2ndDeriv2(nurbsAppr, xi, 'xi');

controlPts = [0,(dy+b_ls/2)/l_ls; 1, 1];
nurbs1D = createNURBSobject(controlPts,[0,0,1,1]);
nurbs1D = elevateNURBSdegree(nurbs1D,1);
nurbs1D = insertKnotsInNURBS(nurbs1D,linspace2(0,1,N(2)-3));
coeffs2(2,1:N(2),2) = nurbs1D.coeffs(1,:);

n = size(coeffs2,2);

p = 2;
t = b_ls/l_ls;
f = @(xi) getNACA(xi,t);
b_u = getNACAapprox(t,p,N,xi_T,f(xi_T));
t = b_us/l_us;
f = @(xi) getNACA(xi,t);
[b_o, Xi] = getNACAapprox(t,p,N,xi_T,f(xi_T));


controlPts = ones(4,n,6);
controlPts(1:2,:,1) = coeffs2(1:2,:,2)*l_ls;
controlPts(3,:,1) = 0;
controlPts(1,:,2) = b_u(:,1)*l_ls;
controlPts(2,:,2) = b_u(:,2)*l_ls;
controlPts(3,:,2) = 0;
controlPts(1,:,3) = delta_s + b_o(:,1)*l_us;
controlPts(2,:,3) = b_o(:,2)*l_us;
controlPts(3,:,3) = h_s;

controlPts(1,:,:) = -controlPts(1,:,:);

controlPts(:,:,4:6) = controlPts(:,:,3:-1:1);
controlPts(2,:,4:6) = -controlPts(2,:,4:6);
controlPts(1,:,6) = controlPts(1,:,1);
controlPts(2,N(2):end,6) = -(dy2+b_ls/2);
controlPts(2,1:N(2),6) = -controlPts(2,1:N(2),1)*(dy2+b_ls/2)/(dy+b_ls/2);


% avgKnot1 = (sqrt(0.3)-xi_t)/(1-xi_t)*aveknt([0,0,0,1,2,3,3,3]/3,3);
% avgKnot2 = (sqrt(0.3)-xi_t)/(1-xi_t)+(1-(sqrt(0.3)-xi_t)/(1-xi_t))*aveknt([0,0,0,1,2,3,4,5,6,6,6]/6,3);
% T = [avgKnot1(2:end-1),avgKnot2];
% P1 = controlPts(1:3,N(2),end);
% P2 = controlPts(1:3,end,end);
% controlPts(1:3,N(2)+1:end,end) = P1*ones(size(T))+(P2-P1)*T;


controlPts = transformPts(controlPts,alpha,x_0);
Eta = [0,0,1,2,3,4,5,5]/5;

nurbs = createNURBSobject(controlPts,{Xi,Eta});

% plotControlPts
% for i = 2:N(2)-1
%     xi = xi_arr(i);
%     P_C = coeffs2(1:2,i,1);
%     g = @(xi) dot(dP(xi),P(xi)-P_C);
%     dgdxi = @(xi) dot(d2P(xi),P(xi)-P_C)+dot(dP(xi),dP(xi));
%     xi = newtonsMethod(g,dgdxi,xi,maxIter,Eps);
%     normal = normalVec(xi,dP);
%     P_T = P(xi);
%     s = (-dx/l_ls - P_T(1))/normal(1);
% %     v = normal*s + P_T;
%     coeffs2(2,i,2) = normal(2)*s + P_T(2);
% end
% temp = coeffs2(1:3,:,2);
% temp(2,N(2):end) = (dy2+b_ls/2)/l_ls;
% temp(1,N(2)+1:N(3)-2) = coeffs2(1,N(2)+1:N(3)-2,1);
% for i = N(2)+1:N(3)-2
%     xi = xi_arr(i);
%     P_C = coeffs2(1:2,i,1);
%     g = @(xi) dot(dP(xi),P(xi)-P_C);
%     dgdxi = @(xi) dot(d2P(xi),P(xi)-P_C)+dot(dP(xi),dP(xi));
%     xi = newtonsMethod(g,dgdxi,xi,maxIter,Eps);
%     P_T = P(xi);
%     normal = normalVec(xi,dP);
%     s = ((dy+b_ls/2)/l_ls - P_T(2))/normal(2);
% %     v = normal*s + P_T;
%     coeffs2(1,i,2) = normal(1)*s + P_T(1);
%     
%     s = ((dy2+b_ls/2)/l_ls - P_T(2))/normal(2);
% %     v = normal*s + P_T;
%     temp(1,i) = normal(1)*s + P_T(1);
% end
% temp = temp*l_ls;
% nurbs3 = createNURBSobject(coeffs2*l_ls,{nurbsAppr2.knots,[0,0,1,1]});

% Eta = [0,0,5,5]/5;
% nurbs = createNURBSobject(controlPts(:,:,1:2),{Xi,Eta});
% nurbs = elevateNURBSdegree(nurbs,[0 1]);
function controlPts = transformPts(controlPts,alpha,x_0)

R_x = rotationXaxis(alpha);
ss = size(controlPts);
for i = 1:ss(2)
    for j = 1:ss(3)
        controlPts(1:3,i,j) = x_0.' + R_x*controlPts(1:3,i,j);
    end
end

function R_x = rotationXaxis(alpha)

R_x = [1, 0,           0;
       0, cos(alpha), -sin(alpha);
       0, sin(alpha),  cos(alpha)];

function n = normalVec(xi,dP)
dPdxi = dP(xi);
n = [-dPdxi(2),dPdxi(1)];


