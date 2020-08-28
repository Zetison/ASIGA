function nurbs = getQuarterDiskData(varargin)
% set default values
options = struct('Xi', [0,0,0,1,1,2,2,3,3,4,4,4]/4, ...
                 'R', [1,0],...
                 'parm', 1, ...
                 'phi', 120*pi/180, ...
                 't', 0, ...
                 'r_trans', 0.9, ...
                 'alignWithAxis', 'Zaxis', ...
                 'x_0', [0,0,0]);
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
        nurbs = explodeNURBS(nurbs,2);
        nurbs = nurbs(1);
    case 2
        nurbs = insertKnotsInNURBS(nurbs([1,5]),{{[], 0.5*[1,1], []}, {0.5*[1,1], 0.5*[1,1], []}});
        nurbs = [explodeNURBS(nurbs(1),2),explodeNURBS(explodeNURBS(nurbs(2),2),1)];
        nurbs = nurbs([1,6]);
        nurbs = [flipNURBSparametrization(mirrorNURBS(nurbs(1),'y'),2),nurbs];
        nurbs = rotateNURBS(nurbs,'theta',pi/4);
        nurbs = scaleNURBSweights(nurbs,nurbs{1}.coeffs(end,end,1));
end

alignWithAxis = options.alignWithAxis;
switch alignWithAxis
    case 'Zaxis'
        % do nothing
    case 'Yaxis'
        nurbs = rotateNURBS(nurbs,'theta',-pi/2,'alignWithAxis', 'Zaxis');
        nurbs = rotateNURBS(nurbs,'theta',-pi/2,'alignWithAxis', 'Xaxis');
    case 'Xaxis'
        nurbs = rotateNURBS(nurbs,'theta',pi/2,'alignWithAxis', 'Zaxis');
        nurbs = rotateNURBS(nurbs,'theta',pi/2,'alignWithAxis', 'Yaxis');
end

nurbs = translateNURBS(nurbs,options.x_0);

