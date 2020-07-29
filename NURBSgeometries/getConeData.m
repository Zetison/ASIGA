function nurbs = getConeData(varargin)
options = struct('Xi', [0,0,0,1,1,2,2,3,3,3]/3, ...   % knot vector of revolution parametric direction
                 'R1', 0,...
                 'R2', 1,...
                 'theta',2*pi,...
                 'h', 1, ...
                 'uniformDegree', true);
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end

R1 = options.R1;
R2 = options.R2;
if numel(R1) == 1
    R1 = [0, R1];
end
if numel(R2) == 1
    R2 = [0, R2];
end
h = options.h;
X = [R1(1),0,0; 
     R2(1),0,h;
     R1(2),0,0; 
     R2(2),0,h];
nurbs = getQuadrilateralData('X',X);
nurbs = revolveNURBS(nurbs,options);
if options.uniformDegree
    nurbs = makeUniformNURBSDegree(nurbs);
end

