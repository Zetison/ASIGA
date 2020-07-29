function nurbs = explodedNURBS(varargin)
options = struct('L', 1, ...
                 'x_0', [0,0,0]);     % angle of revolution in degrees
nurbs = varargin{1};
if nargin > 1
    if numel(varargin) > 2
        newOptions = varargin(2:end);
    else
        newOptions = varargin{2};
    end
    options = updateOptions(options,newOptions);
end
nurbs = explodeNURBS(nurbs);

for patch = 1:numel(nurbs)
    d = nurbs{patch}.d;
    x = mean(nurbs{patch}.coeffs(1:d,:),2);
    n = x.' - options.x_0;
    if ~(norm(n) < 10*eps)
        n = n/norm(n);
        nurbs(patch) = translateNURBS(nurbs(patch),options.L*n);
    end
end