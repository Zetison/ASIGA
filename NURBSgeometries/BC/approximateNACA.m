close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find NACA approximation
setBCParameters
b_ls = 2;
l_ls = 12.3;
% b_ls = 2.4;
% l_ls = 13;
% b_ls = 2;
% l_ls = 12.3;
% b_ls = 0.4;
% l_ls = 2.6;
% l_ls = 2.35;
% b_ls = 0.3;
% l_ls = 2.35;
% b_ls = 0.265;
% l_ls = 2.35;
% b_us = 0.22;
t = b_ls/l_ls;
f_old = @(x) 5*t*(a0*sqrt(x)-a1*x-a2*x.^2+a3*x.^3-a4*x.^4);
f = @(x) getNACA(sqrt(x),t);
100*sqrt(integral(@(x)(f_old(x)-f(x)).^2,0,1)/integral(@(x)f_old(x).^2,0,1))
p = 8;

s_t = 0.0286;
st = 0.2;
Xi = [zeros(1,p+1),s_t*ones(1,p),st*ones(1,p),0.3*ones(1,p),linspace2(0.3,1,8),ones(1,p+1)];
xValues = [];
xIndices = [];
idx = 2;
% xValues = [0.2,0.3,0.4].';
% xIndices = [8,9,10];
cntrlPts = getNACAapprox3(t,p,Xi,idx,xValues,xIndices,st);
coeffs = ones(3,size(cntrlPts,1));
coeffs(1:2,:) = cntrlPts.';

nurbsExact = createNURBSobject(coeffs,Xi);

C = @(xi) evaluateNURBS(nurbsExact, xi);
xi_arr = linspace(0,1,1000);
C_arr = zeros(numel(xi_arr),2);
for i = 1:numel(xi_arr)
    C_arr(i,:) = C(xi_arr(i));
end
figure(1)
plot(l_ls*C_arr(:,1),l_ls*C_arr(:,2), l_ls*xi_arr.^2, l_ls*f(xi_arr.^2), l_ls*xi_arr.^2, l_ls*f_old(xi_arr.^2),l_ls*coeffs(1,:),l_ls*coeffs(2,:),'-o')
hold on
legend('B-spline','Exact','Old exact','Control polygon approx')
return
r = 1.1019*t^2;
theta = linspace(pi/2,pi,1000);
plot(l_ls*r*(1+cos(theta)),l_ls*r*sin(theta))
% return
xi_T = 0.168967783470083;
C_T = C(xi_T);
nurbs2 = insertKnotsInNURBS(nurbsExact,xi_T*ones(1,p));
coeffs2 = nurbs2.coeffs;
coeffs2(:,:,2) = nurbs2.coeffs;
dx = 0.2/13;
coeffs2(1,1:9,2) = -dx;
coeffs2(2,9:end,2) = b_ls/2/l_ls;
nurbs3 = createNURBSobject(coeffs2*l_ls,{nurbs2.knots,[0,0,1,1]});
% figure(42)
% plotNURBS(nurbs3,[100 2], 1, getColor(1));
% hold on
% axis equal
% plotControlPts(nurbs3)

% return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = [1, 5, 9, 16];
p = 2;
xi_T = 0.168967783470083;
% x = [0,0.1,0.3,1];
x = [0,xi_T.^2,0.3,1];
xi = [0,xi_T,sqrt(0.3),1];
[cntrlPts, Xi] = getNACAapprox(t,p,N,xi,f(x));

coeffs = ones(3,size(cntrlPts,1));
coeffs(1:2,:) = cntrlPts.';
%%%%%%%%%%%%
% Xi = [zeros(1,p+1), x(6+3*(N-1)), x(6+3*(N-1):end) ones(1,p+1)];
% coeffs = reshape(coeffs,3,numel(coeffs)/3);
nurbsAppr = createNURBSobject(coeffs,Xi);
P = @(xi) evaluateNURBS(nurbsAppr, xi);
xi_arr = linspace(0,1,1000);
P_arr = zeros(numel(xi_arr),2);
for i = 1:numel(xi_arr)
    P_arr(i,:) = P(xi_arr(i));
end
figure(1)
hold on
% axis equal
plot(l_ls*P_arr(:,1),l_ls*P_arr(:,2),l_ls*coeffs(1,:),l_ls*coeffs(2,:),'-o')
legend('B-spline','Exact','Old exact','Sphere','B-spline approx','Control polygon approx')
100*sqrt(integral(@(x)(f_approx(sqrt(x),P)-f(x)).^2,0,1,'ArrayValued',true)/integral(@(x)f(x).^2,0,1,'ArrayValued',true))
return
nurbsAppr2 = createNURBSobject(coeffs,Xi);
coeffs2 = nurbsAppr2.coeffs;
coeffs2(:,:,2) = nurbsAppr2.coeffs;
dx = 0.2/13;
coeffs2(1,1:N(2),2) = -dx;
xi_arr = aveknt(Xi, p+1);
dP = @(xi) evaluateNURBS_deriv2(nurbsAppr, xi, 'xi');
options = optimset('TolFun',1e-15,'TolX',1e-15,'MaxFunEvals',1000,'MaxIter',1000);
for i = 2:N(2)-1
    xi = xi_arr(i);
    P_C = coeffs2(1:2,i,1);
    xi = fminsearchbnd(@(xi) objFun2(xi,@(xi)P(xi),dP,P_C),xi,xi_arr(i-1),xi_arr(i+1),options);
    normal = normalVec(xi,dP);
    P_T = P(xi);
    s = (-dx - P_T(1))/normal(1);
%     v = normal*s + P_T;
    coeffs2(2,i,2) = normal(2)*s + P_T(2);
end
for i = N(2)+1:N(3)-2
    xi = xi_arr(i);
    P_C = coeffs2(1:2,i,1);
    xi = fminsearchbnd(@(xi) objFun2(xi,@(xi)P(xi),dP,P_C),xi,xi_arr(i-1),xi_arr(i+1),options);
    P_T = P(xi);
    normal = normalVec(xi,dP);
    s = (b_ls/2/l_ls - P_T(2))/normal(2);
%     v = normal*s + P_T;
    coeffs2(1,i,2) = normal(1)*s + P_T(1);
end
coeffs2(2,N(2):end,2) = b_ls/2/l_ls;
nurbs3 = createNURBSobject(coeffs2*l_ls,{nurbsAppr2.knots,[0,0,1,1]});
figure(43)
plotNURBS(nurbs3,[100 2], 1, getColor(1));
% plotControlPts(nurbs3)
hold on
axis equal
return


% return
% temp1 = coeffs(:,3:1+N1);
% temp2 = coeffs(:,4+N1:end-1);
% x = [coeffs(2:3,2).', ...
%      temp1(:).', ...
%      coeffs([1,3],N1+3).', ...
%      temp2(:).', ...
%      linspace2(0,xi_T,N1-1), linspace2(xi_T,1,N2-1)];
% noGpts = 9;
% [Q1D, W1D] = gaussQuad(noGpts);

% objFun(x,Q1D,W1D,C,N1,N2,xi_T,C_T)*100
% % return
% options = optimset('TolFun',1e-14,'TolX',1e-14,'MaxFunEvals',1000,'MaxIter',1000,'Display','iter');
% % x = fminsearch(@(x) objFun(x,Q1D,W1D,C,P_T),x,options);
% UB = ones(size(x));
% % UB(end-5) = 0.9;
% x = fminsearchbnd(@(x) objFun(x,Q1D,W1D,C,N1,N2,xi_T,C_T),x,zeros(size(x)),UB,options);

% figure(2)
% if N1 == 1
%     x_temp = 0;
% else
%     x_temp = x(3*N1-3);
% end
% coeffs = [0, 0, 1, ...
%          0, x(1), x(2), ...
%          x(3:3*(N1-1)+2), ...
%          C_T.', 1, ...
%          x(3*N1), C_T(2)+(C_T(2)-x(3*(N1-1)+1))/(C_T(1)-x_temp)*(x(3*N1)-C_T(1)), x(3*N1+1), ...
%          x((3*N1+2):(1+3*(N1+N2-1))), ...
%          1, 0, 1].';
% p = 2;
% Xi_l =  x(2+3*(N1+N2-1):end);
% Xi = [zeros(1,p+1), Xi_l(1:N1-1), xi_T, xi_T, Xi_l(N1:end) ones(1,p+1)];
% coeffs = reshape(coeffs,3,numel(coeffs)/3);
% nurbsAppr2 = createNURBSobject(coeffs,Xi);
% P = @(xi) evaluateNURBS(nurbsAppr2, xi);
% P_arr = zeros(numel(xi_arr),2);
% for i = 1:numel(xi_arr)
%     P_arr(i,:) = P(xi_arr(i));
% end
% plot(l_ls*C_arr(:,1),l_ls*C_arr(:,2), l_ls*P_arr(:,1),l_ls*P_arr(:,2),l_ls*coeffs(1,:),l_ls*coeffs(2,:),'-o')


% 
% function I = objFun(x,Q1D,W1D,C,N1,N2,xi_T,C_T)
% if N1 == 1
%     x_temp = 0;
% else
%     x_temp = x(3*N1-3);
% end
% coeffs = [0, 0, 1, ...
%          0, x(1), x(2), ...
%          x(3:3*(N1-1)+2), ...
%          C_T.', 1, ...
%          x(3*N1), C_T(2)+(C_T(2)-x(3*(N1-1)+1))/(C_T(1)-x_temp)*(x(3*N1)-C_T(1)), x(3*N1+1), ...
%          x((3*N1+2):(1+3*(N1+N2-1))), ...
%          1, 0, 1].';
% coeffs = reshape(coeffs,3,numel(coeffs)/3);
% p = 2;
% Xi_l =  x(2+3*(N1+N2-1):end);
% Xi_l(Xi_l < 0) = 0;
% Xi_l(Xi_l > 1) = 1;
% Xi = [zeros(1,p+1), Xi_l(1:N1-1), xi_T, xi_T, Xi_l(N1:end) ones(1,p+1)];
% nurbs2 = createNURBSobject(coeffs,Xi);
% 
% P = @(xi) evaluateNURBS(nurbs2, xi);
% dPdxi = @(xi) evaluateNURBS_deriv2(nurbs2, xi,'xi');
% I_ref = 0;
% Eps = 1e-14;
% noItrs = 100;
% I = 0;
% xiGe = insertUniform(Xi,2).^2;
% for i_xi = 1:numel(xiGe)-1
%     for i = 1:numel(Q1D)
%         x = parent2ParametricSpace(xiGe(i_xi:i_xi+1), Q1D(i));
%         xi = sqrt(x);
%         X_C = C(xi);
%         y = X_C(2);
%         xit = newtonsMethod(@(xit) f_(xit,P)-x, @(xit) df_(xit,dPdxi),xi,noItrs,Eps,[0,1]);
%         X_P = P(xit);
%         J_2 = 1/2*(xiGe(i_xi+1)-xiGe(i_xi));
%         I = I + (y-X_P(2))^2*J_2*W1D(i);
%         I_ref = I_ref + y^2*J_2*W1D(i);
%     end
% end
% I = sqrt(I/I_ref);
% % I = I/I_ref;
% end
% function I = integrand(x,P,C,dPdxi)
% Eps = 1e-14;
% noItrs = 100;
% xi = sqrt(x);
% X_C = C(xi);
% y = X_C(2);
% xit = newtonsMethod(@(xit) f_(xit,P)-x,@(xit) df_(xit,dPdxi),xi,noItrs,Eps,[0,1]);
% X_P = P(xit);
% I = (y-X_P(2))^2;
% end
%     
% function f = f_(xi,xit,P,C)
% xi(xi < 0) = 0;
% xi(xi > 1) = 1;
% xit(xit < 0) = 0;
% xit(xit > 1) = 1;
% 
% f = dot(C(xi)-P(xit),C(xi)-P(xit));
% end normalVec
function x = f_(xit,P)
xit(xit < 0) = 0;
xit(xit > 1) = 1;
P = P(xit);
x = P(1);
end
function y = f_approx(xi,P)
P = P(xi);
y = P(2);
end

function f = objFun2(xi,P,dP_,P_C)
dP = dP_(xi);
P_T = P(xi);

f = abs(dot(dP,P_T-P_C));
end
function n = normalVec(xi,dP_)
dP = dP_(xi);
n = [-dP(2),dP(1)];
end
% function df = df_(xi,xit,P,C,dPdxi)
% xi(xi < 0) = 0;
% xi(xi > 1) = 1;
% xit(xit < 0) = 0;
% xit(xit > 1) = 1;
% 
% df = -2*dot(dPdxi(xit),C(xi)-P(xit));
% end
function dx = df_(xit,dPdxi)
xit(xit < 0) = 0;
xit(xit > 1) = 1;

dP = dPdxi(xit);
dx = dP(1);
end

% function d2f = d2f_(xi,xit,P,C,dPdxi,d2Pdxi)
% xi(xi < 0) = 0;
% xi(xi > 1) = 1;
% xit(xit < 0) = 0;
% xit(xit > 1) = 1;
% 
% d2f = -2*(dot(d2Pdxi(xit),C(xi)-P(xit)) - dot(dPdxi(xit),dPdxi(xit)));
% end
% 
% figure(50)
% hold on
% xx = linspace(0,1,1000);
% f_arr = zeros(1,numel(xx));
% for i = 1:numel(xx)
%     xi = sqrt(xx(i));
%     X_C = C(xi);
%     y = X_C(2);
%     xit = newtonsMethod(@(xit) f_(xit,P)-xx(i),@(xit) df_(xit,dPdxi),xi,noItrs,Eps,[0,1]);
%     X_P = P(xit);
%     f_arr(i) = (y-X_P(2))^2;
% end
% plot(xx,f_arr)
% 
% xx = linspace(0,1,1000);
% f_arr = zeros(1,numel(xx));
% for i = 1:numel(xx)
%     f_arr(i) = integrand(xx(i),P,C,dPdxi);
% end
% plot(xx,f_arr)

% plot(xx,f_arr)
% % figure(2)
% % xi_arr = linspace(0,1,1000);
% % df_arr = zeros(numel(xi_arr),1);
% % for i = 1:numel(xi_arr)
% %     df_arr(i,:) = df_(xi,xi_arr(i),P,C,dPdxi);
% %     d2f_arr(i,:) = d2f_(xi,xi_arr(i),P,C,dPdxi,d2Pdxi);
% % end
% % plot(xi_arr,df_arr,xi_arr,d2f_arr)
% % hold on