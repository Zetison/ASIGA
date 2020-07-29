function nurbs = getSectorData(varargin)
options = struct('Xi', [0,0,0,1,1,2,2,3,3,3]/3, ...   % knot vector of revolution parametric direction
                 'R', 1,...
                 'theta', 2*pi, ...
                 'uniformDegree', true);     % angle of revolution in degrees
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end

nurbs = getLineData();
nurbs = revolveNURBS(nurbs,options);
if options.uniformDegree
    nurbs = makeUniformNURBSDegree(nurbs);
end

