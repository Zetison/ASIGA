function nurbs = getCylinderData(varargin)
options = struct('Xi', [0,0,0,1,1,2,2,3,3,3]/3, ...   % knot vector of revolution parametric direction
                 'L', 1, ...   % knot vector of revolution parametric direction
                 'R', [1,0],...
                 'd_p', 3, ...
                 'parm', 1, ...
                 'uniformDegree', true);     % angle of revolution in degrees
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
R = options.R;
switch options.d_p
    case 2
        nurbs1 = getArcData(options);
        nurbs1 = ensure3DNURBS(nurbs1);
    case 3
        if numel(R) == 1
            R = [R,0];
            options.R = R;
        end
        nurbs1 = getDiskData(options);
end
nurbs2 = nurbs1;
nurbs2 = translateNURBS(nurbs2,[0,0,options.L]);
nurbs = loftNURBS({nurbs1,nurbs2});
if options.uniformDegree
    nurbs = makeUniformNURBSDegree(nurbs);
end