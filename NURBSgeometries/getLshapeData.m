function nurbs = getLshapeData(varargin)
% set default values
options = struct('a', 1, ...  % knot vector for azimuthal direction
                 'd', 3,...
                 't', 0);
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
a = options.a;
switch options.d
    case 2
        nurbs = getRectangleData();
        nurbs(2) = translateNURBS(nurbs(1),[a,0]);
        nurbs(3) = translateNURBS(nurbs(1),[0,a]);
    case 3
        nurbs = getPrismData();
        nurbs(2) = translateNURBS(nurbs(1),[a,0,0]);
        nurbs(3) = translateNURBS(nurbs(1),[0,a,0]);
        nurbs(4) = translateNURBS(nurbs(1),[0,0,a]);
end



