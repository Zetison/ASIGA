function nurbs = getBeTSSiSmoothM3Data(varargin)
options = struct('R1', 3,...
                 'R2', 5,...
                 't', 1,...
                 'Xi', [0,0,0,1,1,2,2,3,3,4,4,4]/4, ...
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
x = t*(R2-R1)/sqrt(L^2+(R2-R1)^2);
y = t*L/sqrt(L^2+(R2-R1)^2);
theta = asin(x/t);
theta_eta = asin(y/t);
options.theta_eta = theta_eta;
options.Eta = [0,0,0,1,1,1];
options.C = R1 + t;
hSphere1 = rotateNURBS(getEllipsoidData(rmfield(options,'t')),'rotAxis','Zaxis','theta',pi/2);
hSphere1 = translateNURBS(rotateNURBS(hSphere1,'rotAxis','Yaxis','theta',pi/2),[-R1*tan(theta),0,0]);
hSphere1 = mirrorNURBS(hSphere1,'x');
hSphere1 = flipNURBSparametrization(hSphere1,2);

options.C = R2 + t;
options.theta_eta = pi - theta_eta;
hSphere2 = rotateNURBS(getEllipsoidData(rmfield(options,'t')),'rotAxis','Zaxis','theta',pi/2);
hSphere2 = translateNURBS(rotateNURBS(hSphere2,'rotAxis','Yaxis','theta',pi/2),[-L-R2*tan(theta),0,0]);
cone = loftNURBS({subNURBS(hSphere2,'at',[0,0;0,1]),subNURBS(hSphere1,'at',[0,0;1,0])});
nurbs = [hSphere2, cone, hSphere1];

nurbs = makeUniformNURBSDegree(nurbs);
nurbs = explodeNURBS(nurbs);