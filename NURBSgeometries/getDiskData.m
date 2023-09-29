function nurbs = getDiskData(varargin)
% set default values
options = struct('Xi', [0,0,0,1,1,2,2,3,3,3]/3, ...  % knot vector for azimuthal direction
                 'R', [1,0],...
                 'parm', 1, ...
                 't', 0, ...
                 'phi', 120*pi/180, ...
                 's_trans', 0.9, ...
                 'uniformDegree', true);
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
R = options.R;
R1 = options.R(1);
if numel(R) == 1
    R = 0;
else
    R = options.R(2);
end
parm = options.parm;
if R ~= 0 && parm ~=1
    error('This combination does not exist')
end
t = options.t;
% if t == 0
%     d = 2;
% else
%     d = 3;
% end
w = 1/sqrt(2);
switch parm
    case 1
        nurbs1D = getLineData('x',[R,R1]);
        nurbs = revolveNURBS(nurbs1D,options);
    case 2
        r = R1*options.s_trans;
        phi = options.phi;
        x = r*tan(phi/2)/(1+tan(phi/2));
        w2 = w;
        w3 = w2*w2;
        Xi = [0,0,0,1,1,1];
        Eta = [0,0,1,1];
        coeffs = zeros(3,3,2);
        coeffs(:,:,1) = [	 r	 x	 0;
                             0	 x	 r;
                             1   w2	 1];
        coeffs(:,:,2) = [	 R1	 R1	 0;
                             0	 R1	 R1;
                             1   w	 1];
        nurbs = cell(1,5);
        nurbs(1) = createNURBSobject(coeffs,{Xi, Eta});
        nurbs(1) = permuteNURBS(nurbs(1),[2,1]); % ensure upwards pointing normal vector
        nurbs(2) = rotateNURBS(nurbs(1),'theta',pi/2);
        nurbs(3) = rotateNURBS(nurbs(2),'theta',pi/2);
        nurbs(4) = rotateNURBS(nurbs(3),'theta',pi/2);

        coeffs = zeros(3,3,3);
        coeffs(:,:,1) = [	-r	-x	 0;
                             0  -x	-r;
                             1   w2	 1];
        coeffs(:,:,2) = [	-x	 0	 x;
                             x	 0	-x;
                             w2	 w3	 w2];
        coeffs(:,:,3) = [	 0	 x	 r;
                             r	 x	 0;
                             1   w2	 1];
        Xi = [0,0,0,1,1,1];
        Eta = [0,0,0,1,1,1];
        nurbs(5) = createNURBSobject(coeffs,{Xi, Eta});
    case 3
        Xi = [0 0 0 1 1 1];
        Eta = [0 0 0 1 1 1];

        coeffs = zeros(3,3,3);
        w2 = sqrt(2)-1;
        coeffs(:,:,1) = [	 R1	 R1	 0;
                             0  -R1	-R1;
                             1   w	 1];
        coeffs(:,:,2) = [	 R1	 0	-R1;
                             R1	 0	-R1;
                             w   w2	 w];
        coeffs(:,:,3) = [	 0	-R1	-R1;
                             R1	 R1	 0;
                             1   w	 1];

        nurbs = createNURBSobject(coeffs,{Xi, Eta});
        nurbs = permuteNURBS(nurbs,[2,1]); % ensure upwards pointing normal vector
end
if t == 0
    nurbs = ensure2DNURBS(nurbs);
else
    nurbs = extrudeNURBS(nurbs,'extrudeDir',[0,0,t],'flip',t < 0);
end
if options.uniformDegree
    nurbs = makeUniformNURBSDegree(nurbs);
end



