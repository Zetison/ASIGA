function nurbs = getOctSphereData(varargin)
% set default values
options = struct('R', 1, ...
                 'x_0',[0, 0, 0],...
                 'parm', 1, ...
                 't', 0, ...
                 'scaleWeights', true, ...
                 'prec', 'double');
             
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end

options.C = options.R;
switch options.parm
    case 1
        nurbs = getHalfSphereData(options);
        nurbs = nurbs(1);
    case 2
        nurbs = cell(1,3);
        temp = getQuadSphereData(options);
        nurbs(1) = temp(3);
        nurbs(2) = rotateNURBS(nurbs(1),'rotAxis',[1,1,1],'theta',2*pi/3);
        nurbs(3) = rotateNURBS(nurbs(2),'rotAxis',[1,1,1],'theta',2*pi/3);
end
