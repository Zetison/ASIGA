function nurbs = getLineData(varargin)
options = struct('x', 1);

if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
x = options.x;
if numel(x) == 1
    x = [0,x];
end
controlPts = [x; 
              1, 1];
nurbs = createNURBSobject(controlPts,{[0 0 1 1]});