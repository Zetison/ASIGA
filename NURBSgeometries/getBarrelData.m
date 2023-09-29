function nurbs = getBarrelData(varargin)

options = struct('Xi', [0,0,0,1,1,2,2,3,3,3]/3, ...  % knot vector for azimuthal direction
                 'R', 1,...
                 'parm', 1, ...
                 'L', pi, ...
                 't', 0.1, ...
                 'pmlFill', false, ...
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
    nurbsDisk = permuteNURBS(getDiskData('Xi',Xi,'R',R,'parm',parm,'t',0),[2,1]);
    nurbsTop = translateNURBS(nurbsDisk,[0,0,L]);
    nurbsCyl = getCylinderData('Xi',Xi,'R', R, 'parm', parm, 'L', L, 'd_p', d_p);
    nurbs = [nurbsDisk,nurbsCyl,flipNURBSparametrization(nurbsTop,2)];
else
    nurbsDisk = getDiskData('Xi',Xi,'R',R-t,'parm',parm,'t',t);
                 
    if options.pmlFill
        torus = explodeNURBS(getTorusData('Xi', [0 0 0 1 1 2 2 3 3 4 4 4]/4, 'Eta', Xi, 'R', R-t, 'r', t), 1);
        nurbsDisk = [nurbsDisk, torus(2)];
    else
        nurbsDisk = [nurbsDisk, getDiskData('Xi',Xi,'R',[R,R-t],'parm',1,'t',t)];
    end
    nurbsBot = translateNURBS(mirrorNURBS(nurbsDisk,'z'),[0,0,t]);
    nurbsBot = permuteNURBS(nurbsBot,[2,1,3]);
    nurbsTop = translateNURBS(nurbsDisk,[0,0,L-t]);
    if options.pmlFill
        nurbsTop(2) = permuteNURBS(nurbsTop(2),[2,1,3]);
        nurbsTop(2) = flipNURBSparametrization(nurbsTop(2),2);
    else
        nurbsTop(2) = permuteNURBS(nurbsTop(2),[2,3,1]);
    end
    nurbsTop(1) = permuteNURBS(nurbsTop(1),[2,1,3]);
    nurbsTop(1) = flipNURBSparametrization(nurbsTop(1),2);
    nurbsCyl = getCylinderData('Xi',Xi,'R',[R,R-t],'parm',1,'L',L-2*t);
    nurbsCyl = translateNURBS(nurbsCyl,[0,0,t]);
    nurbsCyl = permuteNURBS(nurbsCyl,[2,3,1]);
    nurbs = [nurbsBot,nurbsCyl,nurbsTop];
end
nurbs = translateNURBS(nurbs,[0,0,-L/2]);
nurbs = rotateNURBS(nurbs,'rotAxis','Yaxis');
nurbs = rotateNURBS(nurbs,'rotAxis','Xaxis', 'theta', pi);
if parm == 2
    nurbs = explodeNURBS(nurbs);
end
