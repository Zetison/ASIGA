function nurbs = getHalfDiskData(varargin)
% set default values
options = struct('R', 1, ...
                 'Xi', [0 0 0 1 1 2 2 3 3 4 4 4]/4,...
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

nurbs = getDiskData(options);
switch options.parm
    case 1
        nurbs = explodeNURBS(nurbs);
        nurbs = nurbs([1,4]);
    case 2
        nurbs([1,3,5]) = insertKnotsInNURBS(nurbs([1,3,5]),{{[], 0.5*[1,1], []}, {[], 0.5*[1,1], []}, {0.5*[1,1], [], []}});
        nurbs = [explodeNURBS(nurbs(1),2), nurbs(2), explodeNURBS(nurbs(3),2), explodeNURBS(nurbs(5),1)];
        nurbs = nurbs([2,3,4,6]);
        nurbs = rotateNURBS(nurbs,'theta',-pi/4-pi/2);
        nurbs = scaleNURBSweights(nurbs,nurbs{1}.coeffs(end,end,1));
end
