function nurbs = getBeTSSiM1Data(varargin)
options = struct('R', 3,...
                 'L', 40,...
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
t = options.t;
L = options.L;
temp = getHalfDiskData('R',R,'parm',2);
nurbs = rotateNURBS(getQuadSphereData('R',R,'parm',2),'rotAxis',[1,0,0],'theta',-pi/2);
nurbs(5:8) = rotateNURBS(temp,'rotAxis',[1,0,0],'theta',pi/2);
arc = subNURBS(nurbs([1,3,4,5,7,8]),'at',{[1,0; 0 0], [0,0; 0 1], [0,0; 1 0], [0,0; 1 0], [0,0; 0 1], [0,1; 0 0]});
arc([3,4,6]) = flipNURBSparametrization(arc([3,4,6]));
arc2 = translateNURBS(arc,[-L,0,0]);
nurbs(9:14) = loftNURBS({arc,arc2});
nurbs(15:18) = translateNURBS(rotateNURBS(rotateNURBS(temp,'rotAxis',[0,0,1],'theta',pi/2),'rotAxis',[0,1,0],'theta',-pi/2),[-L,0,0]);

nurbs = translateNURBS(nurbs,[(L-R)/2, -R/2, 0]);