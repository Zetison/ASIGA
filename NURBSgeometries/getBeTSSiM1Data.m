function nurbs = getBeTSSiM1Data(varargin)
options = struct('R', 3,...
                 'L', 40,...
                 'x_0', [0,0,0], ...
                 'parm', 2, ...
                 't', 0.02); 
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
R = options.R;
L = options.L;
parm = options.parm;

switch parm
    case 1
        halfDisk = getHalfDiskData('R',R,'parm',parm);
        halfDisk1 = rotateNURBS(halfDisk,'rotAxis',[1,0,0],'theta',pi/2);
        halfDisk2 = permuteNURBS(translateNURBS(rotateNURBS(halfDisk1,'rotAxis',[0,0,1],'theta',pi/2),[-L,0,0]),[2,1]);
        halfSphere = rotateNURBS(getQuadSphereData('R',R,'parm',parm),'rotAxis',[0,1,0],'theta',pi/2);
        
        arc2 = subNURBS(halfDisk2,'at',{[0,1; 0 1], [1,0; 0 1]});
        arc2 = arc2([1,3,4,2]);
        arc1 = subNURBS([halfDisk1, halfSphere],'at',{[0,0; 0 1], [0,0; 1 0], [0,0; 1 0], [0,0; 1 0]});
        halfCyl = loftNURBS({arc2,arc1});
        halfCyl(1) = flipNURBSparametrization(halfCyl(1),1);
        nurbs = [halfDisk2, halfCyl, halfSphere, halfDisk1];
    case 2
        temp = getHalfDiskData('R',R,'parm',parm);
        nurbs = rotateNURBS(getQuadSphereData('R',R,'parm',parm),'rotAxis',[1,0,0],'theta',-pi/2);
        nurbs(5:8) = rotateNURBS(temp,'rotAxis',[1,0,0],'theta',pi/2);
        arc = subNURBS(nurbs([1,3,4,5,7,8]),'at',{[1,0; 0 0], [0,0; 0 1], [0,0; 1 0], [0,0; 1 0], [0,0; 0 1], [0,1; 0 0]});
        arc([3,4,6]) = flipNURBSparametrization(arc([3,4,6]));
        arc2 = translateNURBS(arc,[-L,0,0]);
        nurbs(9:14) = loftNURBS({arc,arc2});
        nurbs(15:18) = translateNURBS(rotateNURBS(rotateNURBS(temp,'rotAxis',[0,0,1],'theta',pi/2),'rotAxis',[0,1,0],'theta',-pi/2),[-L,0,0]);
end


nurbs = translateNURBS(nurbs,[(L-R)/2, -R/2, 0]);
nurbs = translateNURBS(nurbs,options.x_0);