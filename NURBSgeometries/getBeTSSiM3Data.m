function nurbs = getBeTSSiM3Data(varargin)
options = struct('R1', 3,...
                 'R2', 5,...
                 'parm', 1, ...
                 't', 0.008, ...
                 'L', 41);
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
t = options.t;
L = options.L;
parm = options.parm;
hSphere1 = getHalfSphereData('R',options.R1,'parm',parm,'t',t);
hSphere1 = rotateNURBS(hSphere1,'rotAxis','Yaxis','theta',pi/2);
hSphere2 = getHalfSphereData('R',options.R2,'parm',parm,'t',t);
hSphere2 = rotateNURBS(hSphere2,'rotAxis','Yaxis','theta',-pi/2);
hSphere2 = translateNURBS(hSphere2,[-L,0,0]);
R1 = options.R1-[t,0];
R2 = options.R2-[t,0];
cone = rotateNURBS(getConeData('R1',R1,'R2',R2,'Xi',[0,0,0,1,1,2,2,3,3,4,4,4]/4,'h',L),'rotAxis','Yaxis','theta',-pi/2);
cone = permuteNURBS(cone,[3,1,2]);
nurbs = [hSphere2, cone, hSphere1];
if parm == 2
    nurbs([1:5,7:end]) = rotateNURBS(nurbs([1:5,7:end]),'rotAxis','Xaxis','theta',pi/4);
    nurbs(6) = elevateNURBSdegree(nurbs(6),[0,0,2]);
end
nurbs = explodeNURBS(nurbs);


