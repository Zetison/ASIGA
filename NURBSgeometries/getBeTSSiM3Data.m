function nurbs = getBeTSSiM3Data(varargin)
options = struct('R1', 3,...
                 'R2', 5,...
                 'parm', 1, ...
                 't', 0.008, ...
                 'Xi', [0,0,0,1,1,2,2,3,3,4,4,4]/4, ...
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
if parm == 1
    Xi = options.Xi;
else
	Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
end
hSphere1 = getHalfSphereData('R',options.R1,'parm',parm,'t',t,'Xi',Xi);
hSphere1 = rotateNURBS(rotateNURBS(hSphere1,'rotAxis','Zaxis','theta',pi),'rotAxis','Yaxis','theta',pi/2);
hSphere2 = getHalfSphereData('R',options.R2,'parm',parm,'t',t,'Xi',Xi);
hSphere2 = mirrorNURBS(hSphere2,'z');
hSphere2 = rotateNURBS(rotateNURBS(hSphere2,'rotAxis','Zaxis','theta',pi),'rotAxis','Yaxis','theta',pi/2);
hSphere2 = translateNURBS(hSphere2,[-L,0,0]);
hSphere2 = flipNURBSparametrization(hSphere2,2);
R1 = options.R1-[t,0];
R2 = options.R2-[t,0];
cone = rotateNURBS(getConeData('R1',R1,'R2',R2,'Xi',Xi,'h',L),'rotAxis','Yaxis','theta',-pi/2);
cone = permuteNURBS(cone,[3,1,2]);
cone = flipNURBSparametrization(cone,[1,2]);
nurbs = [hSphere2, cone, hSphere1];
if parm == 2
    nurbs([1:5,7:end]) = rotateNURBS(nurbs([1:5,7:end]),'rotAxis','Xaxis','theta',pi/4);
    nurbs(6) = elevateNURBSdegree(nurbs(6),[0,0,2]);
end
nurbs = rotateNURBS(nurbs,'rotAxis','Xaxis','theta',-pi/2);
nurbs = explodeNURBS(nurbs);


