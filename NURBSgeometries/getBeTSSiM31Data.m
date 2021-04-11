function nurbs = getBeTSSiM31Data(varargin)
options = struct('R', 3,...
                 'L2', 40,...
                 'R1', 3,...
                 'R2', 5,...
                 'L', 41,...
                 'x_0', [0,0,0], ...
                 't', 0.008, ...
                 't2', 0.02); 
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
L2 = options.L2;
eta = 0.21;
eta2 = L2/L;
shift = L-L2;
halfCyl = permuteNURBS(rotateNURBS(rotateNURBS(getHalfDiskData('R',R,'parm',1,'t',L2),'rotAxis',[1,0,0],'theta',-pi/2),'rotAxis',[0,0,1],'theta',pi/2),[2,3,1]);
halfCyl = flipNURBSparametrization(makeUniformNURBSDegree(halfCyl),[1,2]);
quadSphere = rotateNURBS(getQuadSphereData('R',R,'parm',1,'t',R),'rotAxis',[0,1,0],'theta',pi/2);
quadSphere = insertKnotsInNURBS(quadSphere,{[],eta*[1,1],[]});
M1 = translateNURBS([halfCyl,quadSphere],-[shift,0,0]);
M1rot = rotateNURBS(M1,'rotAxis',[1,0,0],'theta',pi);
halfSphere = flipNURBSparametrization(translateNURBS(rotateNURBS(getHalfSphereData('R',R,'parm',1,'t',R),'rotAxis',[0,1,0],'theta',-pi/2),[-L2-shift,0,0]),[1,2]);
outer = [subNURBS(M1,'at',[0,0; 0 0; 0 1]),subNURBS([halfSphere,M1rot],'at',[0,0; 0 0; 0 1])];
outer = outer([5,8,7,6,2,9,10,1,4,11,12,3]);
outer = explodeNURBS(makeUniformNURBSDegree(outer));
M3options = options;
M3options.x_0 = -[L/2+(options.R2-options.R1)/2,0,0];
M3 = subNURBS(getBeTSSiM3Data(M3options),'at',[0,0; 0 0; 1 0]);
M3(5:8) = insertKnotsInNURBS(M3(5:8),{[],eta2*[1,1],[]});
M3 = explodeNURBS(M3);
M3(5:12) = M3([5,7,9,11,6,8,10,12]);
M1rot = explodeNURBS(M1rot);
outer = outer([1:8,9,11,13,15,10,12,14,16]);
nurbs = [M1rot,halfSphere,loftNURBS({outer,M3})];
nurbs = makeUniformNURBSDegree(nurbs);

nurbs = translateNURBS(nurbs,[L/2+(options.R2-options.R1)/2,0,0]);
nurbs = translateNURBS(nurbs,options.x_0);