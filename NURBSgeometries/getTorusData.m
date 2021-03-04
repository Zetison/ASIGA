function nurbs = getTorusData(varargin)
options = struct('Xi', [0,0,0,1,1,2,2,3,3,3]/3, ...    % knot vector of revolution in the xi-direction
                 'Eta', [0,0,0,1,1,2,2,3,3,3]/3, ...   % knot vector of revolution in the eta-direction
                 'R', 2,...
                 'r', 1,...
                 'parm', 1,...
                 'theta', 2*pi, ...   % angle of revolution in the xi direction
                 'theta_eta', 2*pi);     % angle of revolution in the eta direction 
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
diskOptions = options;
diskOptions.R = options.r;
nurbs = getDiskData(diskOptions);
nurbs = translateNURBS(nurbs,[options.R,0,0]);
nurbs = rotateNURBS(nurbs,'theta',pi/2,'rotAxis','Xaxis');
nurbs = revolveNURBS(nurbs,'Xi',options.Eta,'theta',options.theta_eta);
nurbs = permuteNURBS(nurbs,[2,3,1]);
