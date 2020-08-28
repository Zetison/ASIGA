function nurbs = getQuadSphereData(varargin)
% set default values
options = struct('R', 1, ...
                 'x_0',[0, 0, 0],...
                 'parm', 1, ...
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
        nurbs = getHalfSphereData(options);
        nurbs = nurbs(1:2);
    case 2
        nurbs = cell(1,4);
        temp = getHalfSphereData(options);
        nurbs(2) = temp(2);
        
        temp2 = temp(1);
        temp2 = insertKnotsInNURBS(temp2,{0.5*ones(1,4), [], []});
        temp2 = explodeNURBS(temp2);
        nurbs(1) = temp2(2);
        
        temp2 = temp(3);
        temp2 = insertKnotsInNURBS(temp2,{[], 0.5*ones(1,4), []});
        temp2 = explodeNURBS(temp2);
        nurbs(3) = temp2(1);
        
        temp2 = temp(5);
        temp2 = insertKnotsInNURBS(temp2,{[], 0.5*ones(1,4), []});
        temp2 = explodeNURBS(temp2);
        nurbs(4) = temp2(2);
        nurbs = scaleNURBSweights(nurbs,nurbs{4}.coeffs(end,end,1));
end
