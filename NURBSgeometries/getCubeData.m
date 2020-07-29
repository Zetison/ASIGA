function nurbs = getCubeData(varargin)

options = struct('sideLength', 1, ...
                 'x_0', [0,0,0]);
                   
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end

nurbs = getPrismData('x_0',-[1,1,1]/2);
nurbs = scaleNURBS(nurbs,options.sideLength);
nurbs = translateNURBS(nurbs,options.x_0);