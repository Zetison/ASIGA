function g = sphereMappingFarField(x, y, R, d, type)

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

Y = R*[cos(y).*cos(x), cos(y).*sin(x), sin(y)];

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

if strcmp(type,'g')
    g = Y*d.';
    g = reshape(g,n,m);
    return
end

dYdx = R*[-cos(y).*sin(x), cos(y).*cos(x), zeros(N,1,class(x))];
dYdy = R*[-sin(y).*cos(x),-sin(y).*sin(x), cos(y)];

if strcmp(type,'dgdx')
    g = dYdx*d.';
    g = reshape(g,n,m);
    return
elseif strcmp(type,'dgdy')
    g = dYdy*d.';
    g = reshape(g,n,m);
    return
end

dYdx2 = [-Y(:,1:2),  zeros(N,1,class(x))];
dYdy2 = -Y;
switch type
    case 'dgdx2'
        g = dYdx2*d.';
        g = reshape(g,n,m);
        return
    case 'dgdy2'
        g = dYdy2*d.';
        g = reshape(g,n,m);
        return
    case 'Hessian'
        dYdxdy = R*[sin(y).*sin(x), -sin(y).*cos(x), zeros(N,1,class(x))];
        dgdx2 = dYdx2*d.';
        dgdy2 = dYdy2*d.';
        dgdxdy = dYdxdy*d.';
        g = dgdx2.*dgdy2 - dgdxdy.^2;
        g = reshape(g,n,m);
        return
end
dYdx3 = -dYdx;
dYdy3 = -dYdy;
switch type
    case 'dgdx3'
        g = dYdx3*d.';
        g = reshape(g,n,m);
        return
    case 'dgdy3'
        g = dYdy3*d.';
        g = reshape(g,n,m);
        return
end
dYdx4 = -dYdx2;
dYdy4 = -dYdy2;
switch type
    case 'dgdx4'
        g = dYdx4*d.';
        g = reshape(g,n,m);
        return
    case 'dgdy4'
        g = dYdy4*d.';
        g = reshape(g,n,m);
        return
end
