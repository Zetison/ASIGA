function nurbsCol = addSailToBeTSSi_exact(a0, a1, a2, a3, a4, l_ls, b_ls, l_us, b_us, h_s, delta_s, x_0, alpha, filename, plotNURBS)

t = b_ls/l_ls;
p = 8;
Xi = [zeros(1,p+1) ones(1,p+1)];
n = length(Xi)-(p+1);

f = @(x) 5*t*(a0*sqrt(x)-a1*x-a2*x.^2+a3*x.^3-a4*x.^4);

xi_arr = linspace(0,1,p+1);
A = zeros(p+1);
F = zeros(p+1,2);
for i = 1:p+1
    i1 = findKnotSpan(n, p, xi_arr(i), Xi);
    A(i,:) = Bspline_basis(i1, xi_arr(i), p, Xi, 0);
    x = xi_arr(i)^2;
    F(i,:) = [x, f(x)];
end
b_u = A\F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = b_us/l_us;
f = @(x) 5*t*(a0*sqrt(x)-a1*x-a2*x.^2+a3*x.^3-a4*x.^4);

xi_arr = linspace(0,1,p+1);
A = zeros(p+1);
F = zeros(p+1,2);
for i = 1:p+1
    i1 = findKnotSpan(n, p, xi_arr(i), Xi);
    A(i,:) = Bspline_basis(i1, xi_arr(i), p, Xi, 0);
    x = xi_arr(i)^2;
    F(i,:) = [x, f(x)];
end
b_o = A\F;
% controlPts = [b.'; ones(1,p+1)];

controlPts = ones(4,p+1,2);
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
nurbsCol = {nurbs};
if nargin == 14
    printCtrlPtsToFile(nurbs.coeffs, [filename, '1.txt'])
end
a = 7;
t = b_ls/l_ls;
f = @(x) 5*t*(a0*sqrt(x)-a1*x-a2*x.^2+a3*x.^3-a4*x.^4);
varCol.colorFun = @(X) log10((abs(abs(X(2))/l_ls - f(-(X(1) + 19 - a)/l_ls)))/1.2);

if nargin == 14
    plotNURBS(nurbs,[100 10], 1, 1.5*[44 77 32]/255, 1, NaN, varCol);
    % plotNURBS(nurbs,[500 0], 1, 1.5*[44 77 32]/255, 1);
end
controlPts(2,:,:) = -controlPts(2,:,:);
controlPts2 = transformPts(controlPts,alpha,x_0);
nurbs = createNURBSobject(controlPts2,{Xi,Eta});
nurbsCol{2} = nurbs;
if nargin == 14
    printCtrlPtsToFile(nurbs.coeffs, [filename, '4.txt'])
end
if nargin == 14
    hold on
    plotNURBS(nurbs,[100 10], 1, 1.5*[44 77 32]/255, 1, NaN, varCol);
end
controlPts(1,:,1) = controlPts(1,:,2);
controlPts(2,:,1) = 0;
controlPts(3,:,1) = h_s;
controlPts2 = transformPts(controlPts,alpha,x_0);
nurbs = createNURBSobject(controlPts2,{Xi,Eta});
nurbsCol{3} = nurbs;
if nargin == 14
    printCtrlPtsToFile(nurbs.coeffs, [filename, '2.txt'])
end
if nargin == 14
    plotNURBS(nurbs,[500 0], 1, 1.5*[44 77 32]/255, 1);
end
controlPts(2,:,:) = -controlPts(2,:,:);
controlPts2 = transformPts(controlPts,alpha,x_0);
nurbs = createNURBSobject(controlPts2,{Xi,Eta});
nurbsCol{4} = nurbs;
if nargin == 14
    printCtrlPtsToFile(nurbs.coeffs, [filename, '3.txt'])
end
if nargin == 14
    hold on
    plotNURBS(nurbs,[500 0], 1, 1.5*[44 77 32]/255, 1);
end
% nurbs = leastSquares1D(Xi,p_xi,f)

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
