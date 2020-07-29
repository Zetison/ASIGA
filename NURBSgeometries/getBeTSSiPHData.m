function nurbs = getBeTSSiPHData(varargin)
options = struct('Xi', [0,0,0,1,1,2,2,3,3,3]/3, ...   % knot vector of revolution parametric direction
                 'R1', 3.5, ...  
                 'R2', 4, ...  
                 't',  0.04, ...  
                 'gd', 40, ...
                 'uniformDegree', true);     % angle of revolution in degrees
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
nurbs_o = getPHsurface(options);
t = options.t;
options.R1 = options.R1-t;
options.R2 = options.R2-t;
nurbs_i = getPHsurface(options);
nurbs = loftNURBS({nurbs_i,nurbs_o});
if options.uniformDegree
    nurbs = makeUniformNURBSDegree(nurbs);
end

function nurbs = getPHsurface(options)
R1 = options.R1;
R2 = options.R2;
L = options.gd;

nurbs = getArcData('theta',asin(R1/R2),'Xi',[0,0,0,1,1,1],'R',R2);
nurbs = translateNURBS(nurbs,[-sqrt(R2^2-R1^2),0,0]);
nurbs = [translateNURBS(mirrorNURBS(nurbs,'x'), [-L,0,0]), ...
         translateNURBS(getLineData('x',[-L,0]), [0,R1,0]), ...
         flipNURBSparametrization(nurbs,1)];
nurbs = glueNURBS(nurbs,'xi');
nurbs = revolveNURBS(nurbs,'rotAxis',[1,0,0],'Xi',options.Xi);
nurbs = permuteNURBS(nurbs,[2,1]);
