function g = sphereMapping(x, y, R, X_sp, d, type)

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

Y = R*[sin(y).*cos(x), sin(y).*sin(x), cos(y)];

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
% r = norm2(XmY);
r = sqrt(dot(XmY,XmY,2));
if strcmp(type,'g')
    g = r + Y*d.';
    g = reshape(g,n,m);
    return
end

dYdx = R*[-sin(y).*sin(x), sin(y).*cos(x),  zeros(N,1)];
dYdy = R*[ cos(y).*cos(x), cos(y).*sin(x), -sin(y)];
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

dYdx2 = [-Y(:,1:2),  zeros(N,1)];
dYdy2 = -Y;
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
    case 'Hessian'
        dYdxdy = R*[-cos(y).*sin(x), cos(y).*cos(x), zeros(N,1)];
        drdxdy = -drdy.*drdx./r + (dot(dYdx,dYdy,2) - dot(XmY,dYdxdy,2))./r;
        dgdx2 = drdx2 + dYdx2*d.';
        dgdy2 = drdy2 + dYdy2*d.';
        dgdxdy = drdxdy + dYdxdy*d.';
        g = dgdx2.*dgdy2 - dgdxdy.^2;
        g = reshape(g,n,m);
        return
end
dYdx3 = -dYdx;
dYdy3 = -dYdy;
drdx3 = -3*drdx.*drdx2./r + (3*dot(dYdx,dYdx2,2)-dot(XmY,dYdx3,2))./r;
drdy3 = -3*drdy.*drdy2./r + (3*dot(dYdy,dYdy2,2)-dot(XmY,dYdy3,2))./r;
% drdx3 = -3*drdx.*drdx2./r + real(3*dot(dYdx,dYdx2,2)-dot(XmY,dYdx3,2))./r;
% drdy3 = -3*drdy.*drdy2./r + real(3*dot(dYdy,dYdy2,2)-dot(XmY,dYdy3,2))./r;
switch type
    case 'dgdx3'
        g = drdx3 + dYdx3*d.';
        g = reshape(g,n,m);
        return
    case 'dgdy3'
        g = drdy3 + dYdy3*d.';
        g = reshape(g,n,m);
        return
end
dYdx4 = -dYdx2;
dYdy4 = -dYdy2;
drdx4 = -3*drdx2.^2./r - 4*drdx.*drdx3./r ...
          + (3*dot(dYdx2,dYdx2,2)+4*dot(dYdx,dYdx3,2) - dot(XmY,dYdx4,2))./r;
drdy4 = -3*drdy2.^2./r - 4*drdy.*drdy3./r ...
          + (3*dot(dYdy2,dYdy2,2)+4*dot(dYdy,dYdy3,2) - dot(XmY,dYdy4,2))./r;
% drdx4 = -3*drdx2.^2./r - 4*drdx.*drdx3./r ...
%           + real(3*dot(dYdx2,dYdx2,2)+4*dot(dYdx,dYdx3,2) - dot(XmY,dYdx4,2))./r;
% drdy4 = -3*drdy2.^2./r - 4*drdy.*drdy3./r ...
%           + real(3*dot(dYdy2,dYdy2,2)+4*dot(dYdy,dYdy3,2) - dot(XmY,dYdy4,2))./r;
switch type
    case 'dgdx4'
        g = drdx4 + dYdx4*d.';
        g = reshape(g,n,m);
        return
    case 'dgdy4'
        g = drdy4 + dYdy4*d.';
        g = reshape(g,n,m);
        return
end

function z = dot(x,y,dim)
z = sum(x.*y,dim);
