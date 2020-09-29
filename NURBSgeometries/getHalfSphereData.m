function nurbs = getHalfSphereData(varargin)
% set default values
options = struct('R', 1, ...
                 'x_0',[0, 0, 0],...
                 'parm', 1, ...
                 'Xi', [0,0,0,1,1,2,2,3,3,4,4,4]/4, ...
                 't', 0, ...
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
        nurbs = getEllipsoidData(options);
        nurbs = explodeNURBS(nurbs);
        nurbs = nurbs(2:2:end);
    case 2
        nurbs = cell(1,5);
        temp = getEllipsoidData(options);
        nurbs(1) = temp(1);
        temp = insertKnotsInNURBS(temp(2),{0.5*ones(1,4), [], []});
        temp = explodeNURBS(temp(1),1);
        nurbs(2) = temp(1);
        nurbs(3) = rotateNURBS(nurbs(2));
        nurbs(4) = rotateNURBS(nurbs(3));
        nurbs(5) = rotateNURBS(nurbs(4));
        nurbs = scaleNURBSweights(nurbs,nurbs{2}.coeffs(4,end,1));
end
