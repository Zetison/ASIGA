function [nurbs,RR,c_xy,c_z] = getBeTSSiSmoothM3Data(varargin)
options = struct('R1', 3,...
                 'R2', 5,...
                 't', 1,...
                 'Xi', [0,0,0,1,1,2,2,3,3,4,4,4]/4, ...
                 'x_0', [0,0,0], ...
                 'L', 41);
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
t = options.t;
L = options.L;
R1 = options.R1;
R2 = options.R2;
% x = t*(R2-R1)/sqrt(L^2+(R2-R1)^2);
% y = t*L/sqrt(L^2+(R2-R1)^2);
% a = -(R2-R1)/L;
% w = 1/sqrt(2);
% coeffs1 = [x,    R1+t,         R1+t;
%            R1+y, a*(R1+t-x)+R1+y, 0;
%            1,    w,            1];
% arc1 = createNURBSobject(coeffs1,[0,0,0,1,1,1]);
% coeffs2 = [-L-R2-t, -L-R2-t,         -L+x;
%            0,       a*(-L-R2-t-x)+R1+y, R2+y;
%            1,       w,               1]; 
% RR = [R2+y,R1+y];
dx = t*(R2-R1)/sqrt(L^2+(R2-R1)^2);
dy = t*L/sqrt(L^2+(R2-R1)^2);
c = L/2+(R2+R1)/2;
x1 = c-R1;
y1 = R1;
x2 = -c+R2;
y2 = R2;
x3 = x1+dx;
y3 = y1+dy;
x4 = x2+dx;
y4 = y2+dy;
a = -(R2-R1)/L;

c_z = (L+R1+R2)/2+t;
x = a*c_z^2/(a*x3-y3);
y = a*(x-x3) + y3;
c_xy = sqrt(-a*y*c_z^2/x);
% theta = acos(x1/t);
% w1 = cos(theta/2);
w1 = 1/sqrt(2);
coeffs1 = [x3,      c+t,        c+t;
           y3,   a*(c+t-x3)+y3, 0;
           1,    w1,            1];
arc1 = createNURBSobject(coeffs1,[0,0,0,1,1,1]);

% theta = acos(x2/t);
% w2 = cos(theta/2);
w2 = 1/sqrt(2);
coeffs2 = [-c-t,       -c-t,         x4;
           0,       a*(-c-t-x3)+y3,  y4;
           1,       w2,                 1]; 
RR = [R2+dy,R1+dy];
arc2 = createNURBSobject(coeffs2,[0,0,0,1,1,1]);
line = createNURBSobject([coeffs2(:,end), coeffs1(:,1)],[0,0,1,1]);
nurbs = uniteNURBS({arc2,line,arc1});
nurbs = revolveNURBS(nurbs,'rotAxis', [1, 0, 0], 'Xi', options.Xi);
nurbs = permuteNURBS(nurbs,[2,1]);
nurbs = makeUniformNURBSDegree(nurbs);
nurbs = explodeNURBS(nurbs);
% nurbs = translateNURBS(nurbs,[L/2+(R2-R1)/2,0,0]);
nurbs = translateNURBS(nurbs,options.x_0);