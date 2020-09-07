function nurbs = getBarrelData(varargin)

options = struct('Xi', [0,0,0,1,1,2,2,3,3,3]/3, ...  % knot vector for azimuthal direction
                 'R', 1,...
                 'parm', 1, ...
                 'L', pi, ...
                 't', 0.1, ...
                 'd_p', 3, ...
                 'phi', 120*pi/180, ...
                 's_trans', 0.9);
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
parm = options.parm;
R = options.R;
L = options.L;
t = options.t;
d_p = options.d_p;
if parm == 2
    Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4; % p = 2
else
    Xi = options.Xi;
end
if d_p == 2
    nurbsDisk = getDiskData('Xi',Xi,'R',R,'parm',parm,'t',0);
    nurbsTop = translateNURBS(nurbsDisk,[0,0,L]);
    nurbsCyl = getCylinderData('Xi',Xi,'R', R, 'parm', parm, 'L', L, 'd_p', d_p);
    nurbs = [rotateNURBS(nurbsDisk,'rotAxis','Yaxis','theta',pi),nurbsCyl,nurbsTop];
else
    nurbsDisk = [getDiskData('Xi',Xi,'R',R-t,'parm',parm,'t',t), ...
                 getDiskData('Xi',Xi,'R',[R,R-t],'parm',1,'t',t)];
    nurbsTop = translateNURBS(nurbsDisk,[0,0,L-t]);
    nurbsCyl = getCylinderData('Xi',Xi,'R',[R,R-t],'parm',1,'L',L-2*t);
    nurbsCyl = translateNURBS(nurbsCyl,[0,0,t]);
    nurbs = [nurbsDisk,nurbsCyl,nurbsTop];
end
nurbs = explodeNURBS(nurbs);
nurbs = rotateNURBS(nurbs,'rotAxis','Yaxis');
