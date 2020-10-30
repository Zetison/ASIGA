function nurbs = getScordelisLoRoofData(varargin)
options = struct('Xi', [0,0,0,1,1,1], ...   % knot vector of revolution parametric direction
                 'R', 25,...
                 'L', 50, ...   % knot vector of revolution parametric direction
                 't', 0.25, ...
                 'degree', [2,2,2], ...
                 'phi', 40);     % angle of revolution in degrees
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
R = options.R;
t = options.t;
L = options.L;
nurbs = getCylinderData('Xi', options.Xi, 'R', [R+t/2,R-t/2], 'L', L/2, 'theta', options.phi*pi/180);
nurbs = rotateNURBS(nurbs,'rotAxis','Zaxis','theta',-pi/2);
nurbs = rotateNURBS(nurbs,'rotAxis','Xaxis','theta',-pi/2);
nurbs = makeUniformNURBSDegree(nurbs,options.degree);