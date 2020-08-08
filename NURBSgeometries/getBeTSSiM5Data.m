function nurbs = getBeTSSiM5Data(varargin)

options = struct('Xi', [0,0,0,1,1,2,2,3,3,3]/3, ...  % knot vector for azimuthal direction
                 'type','B',...
                 'R', 0.6/2,...
                 'parm', 1, ...
                 'L', 4.8, ...
                 'l', 1, ...
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
l = options.l;
if parm == 2
    Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4; % p = 2
else
    Xi = options.Xi;
end
nurbsDisk = getDiskData('Xi',Xi,'R',R,'parm',parm);
nurbsTop = translateNURBS(nurbsDisk,[0,0,L]);
nurbsCyl = getCylinderData('Xi',Xi,'R',R,'parm',1,'L',L,'d_p',2);
nurbs = [nurbsDisk,nurbsCyl,nurbsTop];
nurbs = [translateNURBS(nurbs,[l/2,0,0]),...
         translateNURBS(nurbs,[-l/2,0,0])];
if strcmp(options.type,'A')
    nurbs = rotateNURBS(nurbs,'rotAxis','Yaxis');
end
