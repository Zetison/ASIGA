function nurbs = getBeTSSiM4Data(varargin)
options = struct('Xi', [0,0,0,1,1,2,2,3,3,4,4,4]/4, ...
                 'parm', 1, ...
                 'phi', 120*pi/180, ...
                 't', 0.02, ...
                 'R', 3, ...
                 'r_trans', 0.9, ...
                 'uniformDegree', true);
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
t = options.t;
R = options.R;
parm = options.parm;
nurbs = getQuarterDiskData('R',R,'alignWithAxis','Zaxis','x_0',[1,1,1]*t/2,'parm',parm,'t',-t,'r_trans',options.r_trans);
if parm == 2
    nurbs(4) = extrudeNURBS(subNURBS(nurbs(2),'at',{[0 0; 0 1; 0 0]}), 'extrudeDir',[-t,0,0], 'flip', false);
    nurbs(5) = extrudeNURBS(subNURBS(nurbs(3),'at',{[1 0; 0 0; 0 0]}), 'extrudeDir',[-t,0,0], 'flip', true);
else
    nurbs(2) = extrudeNURBS(subNURBS(nurbs,'at',{[0 0; 0 1; 0 0]}), 'extrudeDir',[-t,0,0]);
end
n = numel(nurbs);
nurbs(n+1:2*n) = rotateNURBS(nurbs(1:n),'rotAxis',[1,1,1],'theta',2*pi/3);
nurbs(2*n+1:3*n) = rotateNURBS(nurbs(1:n),'rotAxis',[1,1,1],'theta',4*pi/3);
if parm == 1
    nurbs(end+1) = extrudeNURBS(subNURBS(nurbs(n),'at',{[0 0; 1 0; 0 0]}), 'extrudeDir',[0,-t,0], 'flip', true);
else
    nurbs(end+1) = extrudeNURBS(subNURBS(nurbs(n),'at',{[1 0; 0 0; 0 0],[]}), 'extrudeDir',[0,-t,0], 'flip', true);
end
if options.uniformDegree
    nurbs = makeUniformNURBSDegree(nurbs);
end

