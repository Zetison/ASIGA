function g = sphereMappingNURBS(x, y, nurbs, X_sp, d, type)

if numel(x) ~= numel(y)
    if numel(x) == 1
        [n,m] = size(y);
        x = repmat(x,n,m);
    elseif numel(y) == 1
        [n,m] = size(x);
        y = repmat(y,n,m);
    else
        error('The number of elements in x and y should be the same');
    end
end
[n,m] = size(y);
N = n*m;
x = reshape(x,N,1);
y = reshape(y,N,1);

[Y, dYdx, dYdy, dYdx2, dYdy2, dYdxdy] = evaluateNURBS_2ndDeriv(nurbs, [x,y]);

if strcmp(type,'Y1')
    g = Y(:,1);
    g = reshape(g,n,m);
    return
elseif strcmp(type,'Y2')
    g = Y(:,2);
    g = reshape(g,n,m);
    return
elseif strcmp(type,'Y3')
    g = Y(:,3);
    g = reshape(g,n,m);
    return
end

XmY = repmat(X_sp,N,1)-Y;

r = sqrt(dot(XmY,XmY,2));
if strcmp(type,'g')
    g = r + Y*d.';
    g = reshape(g,n,m);
    return
end

drdx = -dot(XmY,dYdx,2)./r;
drdy = -dot(XmY,dYdy,2)./r;

if strcmp(type,'dgdx')
    g = drdx + dYdx*d.';
    g = reshape(g,n,m);
    return
elseif strcmp(type,'dgdy')
    g = drdy + dYdy*d.';
    g = reshape(g,n,m);
    return
end

drdx2 = -drdx.^2./r + (dot(dYdx,dYdx,2) - dot(XmY,dYdx2,2))./r;
drdy2 = -drdy.^2./r + (dot(dYdy,dYdy,2) - dot(XmY,dYdy2,2))./r;
% drdx2 = -drdx.^2./r + real(dot(dYdx,dYdx,2) - dot(XmY,dYdx2,2))./r;
% drdy2 = -drdy.^2./r + real(dot(dYdy,dYdy,2) - dot(XmY,dYdy2,2))./r;
switch type
    case 'dgdx2'
        g = drdx2 + dYdx2*d.';
        g = reshape(g,n,m);
        return
    case 'dgdy2'
        g = drdy2 + dYdy2*d.';
        g = reshape(g,n,m);
        return
    case 'dgdxdy'
        drdxdy = -drdy.*drdx./r + (dot(dYdx,dYdy,2) - dot(XmY,dYdxdy,2))./r;
        g = drdxdy + dYdxdy*d.';
        g = reshape(g,n,m);
        return
    case 'Hessian'
        drdxdy = -drdy.*drdx./r + (dot(dYdx,dYdy,2) - dot(XmY,dYdxdy,2))./r;
        dgdx2 = drdx2 + dYdx2*d.';
        dgdy2 = drdy2 + dYdy2*d.';
        dgdxdy = drdxdy + dYdxdy*d.';
        g = dgdx2.*dgdy2 - dgdxdy.^2;
        g = reshape(g,n,m);
        return
end
% dYdxi3 = -dYdxi;
% dYdeta3 = -dYdeta;
% dYdx3 = dYdxi3*dxidx^3;
% dYdy3 = dYdeta3*detady^3;
% % drdx3 = -3*drdx.*drdx2./r + (3*dot(dYdx,dYdx2,2)-dot(XmY,dYdx3,2))./r;
% % drdy3 = -3*drdy.*drdy2./r + (3*dot(dYdy,dYdy2,2)-dot(XmY,dYdy3,2))./r;
% drdx3 = -3*drdx.*drdx2./r + real(3*dot(dYdx,dYdx2,2)-dot(XmY,dYdx3,2))./r;
% drdy3 = -3*drdy.*drdy2./r + real(3*dot(dYdy,dYdy2,2)-dot(XmY,dYdy3,2))./r;
% switch type
%     case 'dgdx3'
%         g = drdx3 + dYdx3*d.';
%         g = reshape(g,n,m);
%         return
%     case 'dgdy3'
%         g = drdy3 + dYdy3*d.';
%         g = reshape(g,n,m);
%         return
% end
function z = dot(x,y,dim)
z = sum(x.*y,dim);

function [v, dvdxi, dvdeta, d2vdxi2, d2vdeta2, d2vdxideta] = evaluateNURBS_2ndDeriv(nurbs, parm_pt)

d = size(nurbs.coeffs,1) - 1;

p = nurbs.degree(1);
q = nurbs.degree(2);

n = nurbs.number(1);
m = nurbs.number(2);

Xi = nurbs.knots{1};
Eta = nurbs.knots{2};

P = nurbs.coeffs;

i1 = findKnotSpan(n, p, xiE(1), Xi);
i2 = findKnotSpan(m, q, etaE(1), Eta);


[N, dNdxi, d2Ndxi2] = Bspline_basisDers2(i1, parm_pt(:,1), p, Xi);
[M, dMdeta, d2Mdeta2] = Bspline_basisDers2(i2, parm_pt(:,2), q, Eta);

noxi = size(parm_pt,1);
W = zeros(noxi,1);
dWdxi = zeros(noxi,1);
dWdeta = zeros(noxi,1);
d2Wdxi2 = zeros(noxi,1);
d2Wdeta2 = zeros(noxi,1);
d2Wdxideta = zeros(noxi,1);

for k2 = 1:q+1
    A2 = i2 - q + k2 - 1;
    for k1 = 1:p+1
        A1 = i1 - p + k1 - 1;
        weight = P(d+1, A1, A2);

        W           = W             + N(:,k1)     .*M(:,k2)     *weight;
        dWdxi       = dWdxi         + dNdxi(:,k1) .*M(:,k2)     *weight;
        dWdeta      = dWdeta        + N(:,k1)     .*dMdeta(:,k2)*weight;
        d2Wdxi2     = d2Wdxi2       + d2Ndxi2(:,k1).*M(:,k2)     *weight;
        d2Wdeta2    = d2Wdeta2      + N(:,k1)     .*d2Mdeta2(:,k2)*weight;
        d2Wdxideta  = d2Wdxideta	+ dNdxi(:,k1) .*dMdeta(:,k2)*weight;
    end
end

v = zeros(noxi,d);
dvdxi = zeros(noxi,d);
dvdeta = zeros(noxi,d);
d2vdxi2 = zeros(noxi,d);
d2vdeta2 = zeros(noxi,d);
d2vdxideta = zeros(noxi,d);
counter = 1;
for k2 = 1:q+1
    A2 = i2 - q + k2 - 1;
    for k1 = 1:p+1
        A1 = i1 - p + k1 - 1;
        weight = P(d+1, A1, A2);   
        point = P(1:d, A1, A2).';   
        fact = weight./(W.*W);

        NM = N(:,k1).*M(:,k2);
        v = v + NM.*W.*fact*point;

        dvdxi   = dvdxi + (dNdxi(:,k1)  .*M(:,k2).*W - NM.*dWdxi).*fact*point;
        dvdeta  = dvdeta + (dMdeta(:,k2) .*N(:,k1).*W - NM.*dWdeta).*fact*point;

        d2vdxi2   = d2vdxi2 + (d2Ndxi2(:,k1).*W - N(:,k1).*d2Wdxi2 - 2*(dNdxi(:,k1).*W - N(:,k1).*dWdxi).*dWdxi./W).*M(:,k2).*fact*point;
        d2vdeta2   = d2vdeta2 + (d2Mdeta2(:,k2).*W - M(:,k2).*d2Wdeta2 - 2*(dMdeta(:,k2).*W - M(:,k2).*dWdeta).*dWdeta./W).*N(:,k1).*fact*point;
        d2vdxideta   = d2vdxideta + (-dNdxi(:,k1).*dWdeta.*M(:,k2).*W - N(:,k1).*dWdxi.*dMdeta(:,k2).*W - NM.*d2Wdxideta.*W ...
                                        + dNdxi(:,k1).*dMdeta(:,k2).*W.^2 + 2*dWdxi.*dWdeta.*NM)./W.*fact*point;

        counter = counter + 1;
    end
end
