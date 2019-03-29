function addSailToBeTSSi_approx(a0, a1, a2, a3, a4, l_ls, b_ls, l_us, b_us, h_s, delta_s, x_0, alpha, filename)

t_u = b_ls/l_ls;
p = 2;
Xi = [zeros(1,p+1), linspace2(0,1,7), ones(1,p+1)];
n = length(Xi)-(p+1);

f = @(x) 5*t_u*(a0*sqrt(x)-a1*x-a2*x.^2+a3*x.^3-a4*x.^4);
dfdx = @(x) 5*(a0*1/(2*sqrt(x))-a1-2*a2*x+3*a3*x.^2-4*a4*x.^3);
xi_arr = linspace(0,1,n);
A = zeros(n);
F = zeros(n,2);
for i = 1:n
    i1 = findKnotSpan(n, p, xi_arr(i), Xi);
    A(i,i1-p:i1) = Bspline_basis(i1, xi_arr(i), p, Xi, 0);
    x = xi_arr(i)^2;
    F(i,:) = [x, f(x)];
end
b_u = A\F;
figure(2)
controlPts = b_u.';
controlPts(3,:) = 1;
nurbs = createNURBSobject(controlPts,Xi);

xi = linspace(0,1,1000);
X = zeros(2,numel(b_u));
for i = 1:numel(xi)
    X(:,i) = evaluateNURBS(nurbs,xi(i));
end
plot(X(1,:),X(2,:))
hold on
plot(xi,f(xi))
figure(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_o = b_us/l_us;
f = @(x) 5*t_o*(a0*sqrt(x)-a1*x-a2*x.^2+a3*x.^3-a4*x.^4);

xi_arr = linspace(0,1,n);
A = zeros(n);
F = zeros(n,2);
for i = 1:n
    i1 = findKnotSpan(n, p, xi_arr(i), Xi);
    A(i,i1-p:i1) = Bspline_basis(i1, xi_arr(i), p, Xi, 0);
    x = xi_arr(i)^2;
    F(i,:) = [x, f(x)];
end
b_o = A\F;
% controlPts = [b.'; ones(1,p+1)];

controlPts = ones(4,n,2);
controlPts(1,:,1) = b_u(:,1)*l_ls;
controlPts(2,:,1) = b_u(:,2)*l_ls;
controlPts(3,:,1) = 0;
controlPts(1,:,2) = delta_s + b_o(:,1)*l_us;
controlPts(2,:,2) = b_o(:,2)*l_us;
controlPts(3,:,2) = h_s;
controlPts(1,:,:) = -controlPts(1,:,:);
controlPts2 = transformPts(controlPts,alpha,x_0);
Eta = [0,0,1,1];
nurbs = createNURBSobject(controlPts2,{Xi,Eta});
if nargin == 14
    printCtrlPtsToFile(nurbs.coeffs, [filename, '1.txt'])
end
a = 7;
c = 4;
t_u = b_ls/l_ls;
f = @(x,t) 5*t*(a0*sqrt(x)-a1*x-a2*x.^2+a3*x.^3-a4*x.^4);
ss = @(z) (z-c)/h_s;

varCol.colorFun = @(X) log10((abs(abs(X(2))/(ss(X(3))*l_us+(1-ss(X(3)))*l_ls) - f(-(X(1) + 19 - a + ss(X(3))*delta_s)/(ss(X(3))*l_us+(1-ss(X(3)))*l_ls), ss(X(3))*t_o+(1-ss(X(3)))*t_u))/1.2));

plotNURBS(nurbs,[200 0], 1, 1.5*[44 77 32]/255, 1, NaN, varCol);
controlPts(2,:,:) = -controlPts(2,:,:);
controlPts2 = transformPts(controlPts,alpha,x_0);
nurbs = createNURBSobject(controlPts2,{Xi,Eta});
if nargin == 14
    printCtrlPtsToFile(nurbs.coeffs, [filename, '4.txt'])
end
hold on
plotNURBS(nurbs,[200 0], 1, 1.5*[44 77 32]/255, 1, NaN, varCol);

controlPts(1,:,1) = controlPts(1,:,2);
controlPts(2,:,1) = 0;
controlPts(3,:,1) = h_s;
controlPts2 = transformPts(controlPts,alpha,x_0);
nurbs = createNURBSobject(controlPts2,{Xi,Eta});
if nargin == 14
    printCtrlPtsToFile(nurbs.coeffs, [filename, '2.txt'])
end
plotNURBS(nurbs,[200 0], 1, 1.5*[44 77 32]/255, 1);
controlPts(2,:,:) = -controlPts(2,:,:);
controlPts2 = transformPts(controlPts,alpha,x_0);
nurbs = createNURBSobject(controlPts2,{Xi,Eta});
if nargin == 14
    printCtrlPtsToFile(nurbs.coeffs, [filename, '3.txt'])
end
hold on
plotNURBS(nurbs,[200 0], 1, 1.5*[44 77 32]/255, 1);


function controlPts = transformPts(controlPts,alpha,x_0)

R_x = rotationXaxis(alpha);
ss = size(controlPts);
for i = 1:ss(2)
    for j = 1:ss(3)
        controlPts(1:3,i,j) = x_0.' + R_x*controlPts(1:3,i,j);
    end
end

function printCtrlPtsToFile(coeffs, filename)

sz = size(coeffs);
M = zeros(sz(2)*sz(3),3);
for i = 1:sz(2)
    for j = 1:sz(3)
        index = i + (j-1)*sz(2);
        M(index,:) = coeffs(1:3,i,j);
    end
end
dlmwrite(filename,M,'precision',16)

function R_x = rotationXaxis(alpha)

R_x = [1, 0,           0;
       0, cos(alpha), -sin(alpha);
       0, sin(alpha),  cos(alpha)];
