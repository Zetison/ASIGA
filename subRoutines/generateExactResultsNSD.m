clear all

npts = 2;
noSubElements = 2000;
digits(100)
prec = 'mp';
Q = zeros(npts^2,2,prec);
W = zeros(npts^2,1,prec);
% Gauss quadrature along two directions

if strcmp(prec,'mp')
    [ptU, wtU] = gaussQuad_mp(npts,1);
    [ptV, wtV] = gaussQuad_mp(npts,1);
else
    [ptU, wtU] = gaussQuad_sym(npts,1);
    [ptV, wtV] = gaussQuad_sym(npts,1);
end

A = 1;
for j = 1:npts
    for i = 1:npts
        Q(A,:) = [ptU(i), ptV(j)];
        W(A)  = wtU(i)*wtV(j);
        A = A+1;
    end
end

d = -[1,2,3];
nOmegas = 1;
omega_start = 0;
omega_end = 3;
if strcmp(prec,'mp')
    d = mp(d);
    PI = mp('pi');
    R = mp('5');
    P_inc = mp('1');
    omega_arr = 10.^linspace(mp(omega_start),mp(omega_end),nOmegas);
    omega_arr = mp('1000');
elseif strcmp(prec,'sym')
    d = vpa(d);
    PI = vpa('pi');
    R = vpa('5');
    P_inc = vpa('1');
    omega_arr = 10.^linspace(vpa(omega_start),vpa(omega_end),nOmegas);
    omega_arr = vpa('1000');
end
     
d = d/norm(d);
xhat = -d;
dmx = 2*d;   
noElems = 8;
elementList = 1:(noElems-1);
% elementList = 5;
index = [copyVector((1:4).',4,1), copyVector((1:4).',4,2)];
arr = PI*([0, 0.5, 1].'-1/2);
elRangeEta = [arr(1:end-1), arr(2:end)];
arr = (2*[0, 0.25, 0.5, 0.75, 1].'-1)*PI;
elRangeXi = [arr(1:end-1), arr(2:end)];

g =    @(x,y) sphereMappingFarField(x, y, R, dmx, 'g');
%         dgdx = @(x,y) sphereMappingFarField(x, y, R, dmx, 'dgdx');
%         f = @(x,y) dgdx(x,y);
xdY = @(xi,eta) R*(xhat(1)*cos(eta).*cos(xi)+xhat(2)*cos(eta).*sin(xi)+xhat(3)*sin(eta));
f = @(xi,eta) heaviside(xdY(xi,eta)).*xdY(xi,eta).*cos(eta);

I = zeros(nOmegas,noElems,prec);
for i_omega = 1:nOmegas
    parfor e = elementList
        omega = omega_arr(i_omega);
        idXi = index(e+1,1);
        idEta = index(e+1,2);
        xiE = elRangeXi(idXi,:);
        etaE = elRangeEta(idEta,:);
%         noSubElements = 400;
        xiGe = linspace(xiE(1),xiE(2),noSubElements);
        etaGe = linspace(etaE(1),etaE(2),noSubElements);
        tic
        for i_eta = 1:numel(etaGe)-1
            for i_xi = 1:numel(xiGe)-1
                x  = parent2ParametricSpace(xiGe(i_xi:i_xi+1), Q(:,1));
                y  = parent2ParametricSpace(etaGe(i_eta:i_eta+1), Q(:,2));
                J2 = (xiGe(i_xi+1)-xiGe(i_xi))*(etaGe(i_eta+1)-etaGe(i_eta))/4;
                fact = J2*W;
                I(i_omega,e+1) = I(i_omega,e+1) + sum(f(x,y).*exp(1i*omega*g(x,y)).*fact);
            end
        end
        toc
    end
end
% I(:,e)
omega = omega_arr;
P_inc = 1;
R = 5;
sum(I,2)-exactKDT(omega,P_inc,R)/(-1i*k*P_inc*R/(2*PI))
% 
% R = 5;
% P_inc = 1;
% d = -[1,2,3];
% xhat = -d;
% dmx = 2*d;
% g2 =    @(x,y) sphereMappingFarField(x, y, R, dmx, 'g');
% integrand2 = @(y,omega) 1/(1i*double(omega))*(exp(1i*double(omega)*g2(double(xiE(2)),y)) - exp(1i*double(omega)*g2(double(xiE(1)),y)));
% I_Exact = @(omega) integral(@(y)integrand2(y,double(omega)),double(etaE(1)),double(etaE(2)),'AbsTol',1e-14,'RelTol',1e-14);
% 
% I_Exact(omega_arr(1))
% 
% 
% 
function y = heaviside(x)

y = x > 0;
end