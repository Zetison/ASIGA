function nurbs = getMockShellData(varargin)
options = struct('Xi', [0,0,0,1,1,2,2,3,3,3]/3, ...   % knot vector of revolution parametric direction
                 'R', 1,...
                 'parm', 1, ...
                 't', 0.2, ...
                 'L', 10);
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
options.R1 = options.R;
options.R2 = options.R;
nurbs = getBeTSSiM3Data(options);