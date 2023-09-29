function nurbs = getHalfDiskData(varargin)
% set default values
options = struct('R', 1, ...
                 'Xi', [0 0 0 1 1 2 2 3 3 4 4 4]/4,...
                 'x_0',[0, 0, 0],...
                 'parm', 1, ...
                 't', 0, ...
                 'theta3', pi/4, ...
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
    case 3
        theta =  options.theta3;
        R = options.R;
        edges = cell(1,4);
        edges(2) = getArcData('Xi', [0,0,0,1,1,1], 'theta', theta,'R',R);
        edges(1) = flipNURBSparametrization(rotateNURBS(edges(2),'theta',pi-theta),1);
        if theta < pi/4
            XI = [0,0,0,0.5,0.5,1,1,1];
            newKnots = [0.5,0.5];
        else
            XI = [0,0,0,1,1,1];
            newKnots = [];
        end
        edges(3) = translateNURBS(getLineData('x',2*R),[-R,0,0]);
        edges(3) = elevateNURBSdegree(edges(3),1);
        edges(3) = insertKnotsInNURBS(edges(3),{newKnots});        
        
        edges(4) = flipNURBSparametrization(rotateNURBS(getArcData('Xi', XI, 'theta', pi-2*theta),'theta',theta,'R',R),1);
        edges = ensure2DNURBS(edges);
        nurbs = GordonHall(edges);
        nurbs = rotateNURBS(nurbs,'theta', -pi/2);
        nurbs = ensure2DNURBS(nurbs);
end
