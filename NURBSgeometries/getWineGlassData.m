function nurbs = getWineGlassData(varargin)
options = struct('uniformDegree', true);     % angle of revolution in degrees
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end

t1 = 0.001;
t2 = 0.0035; % Thickness of base
k = 0.001;

h1 = 0.005;
h2 = 0.115;
h6 = 0.125;
h7 = 0.15; 
h8 = 0.2;
h9 = 0.22;

r = 0.017/2; % Shaft radius at top and bottom

x_1 = 0.6*r;
y_1 = h1+0.6*r;

x_2 = x_1;
y_2 = h2-0.6*r;


R_0i = 0.073/2; % messured
R_1i = 0.02;
R_2i = 0.002;
R_5i = 0.005;   % bottom of cup
R_8i = 0.074/2; % messured
R_6i = R_8i - 0.0;
R_7i = R_8i + 0.017;

R_0o = R_0i;
R_1o = 0.6*R_0o;
R_4o = 0.2*x_1; %bending radius of shaft
R_6o = R_6i+t1;
R_7o = R_7i+t1;
R_8o = R_8i+t1;
n = 13;
Xi = [0 0 0 0.5 1 1 1.1 1.9 2 2 2.2 2.5 3.2 4 4 4]/4;
controlPts = [R_0i, R_1i,   R_2i,   0, 0,       0,          0,      0, 	R_5i,	R_6i,       R_7i,   R_8i,   R_8i;
              zeros(1,n);
              0,    0,      h1,     h1,y_1+k,   2*(h1+h2)/4,y_2-k,  h2,	h2,     h6+t1/2,	h7,     h8,     h9;
              ones(1,n)];
nurbs_i = createNURBSobject(controlPts,{Xi});
controlPts = [R_0o, R_1o,	x_1+k,	x_1,    x_1-k,	R_4o,   	x_2-k,	x_2, 	x_2+k,	R_6o, 	R_7o,   R_8o,   R_8o;
              zeros(1,n);
              t2,	t2,  	y_1-k,	y_1,    y_1+k,	(h1+h2)/2,  y_2-k,	y_2,	y_2+k,	h6,     h7,     h8,     h9;
              ones(1,n)];
nurbs_o = createNURBSobject(controlPts,{Xi});
nurbs = loftNURBS({nurbs_i,nurbs_o});
nurbs = insertKnotsInNURBS(nurbs,{[] 0.5 []});
nurbs = insertKnotsInNURBS(nurbs,{[0.5, 1.1, 1.9]/4 [] []});

options.Xi = [0,0,0,1,1,2,2,2]/2;
options.theta = pi;
nurbsArcTop = getSectorData(options);
nurbsArcBot = scaleNURBS(nurbsArcTop,t2/2);
nurbsArcBot = rotateNURBS(nurbsArcBot,'theta',pi/2,'rotAxis','Xaxis');
nurbsArcBot = rotateNURBS(nurbsArcBot,'theta',pi/2,'rotAxis','Yaxis');
nurbsArcBot = translateNURBS(nurbsArcBot,[R_0o,0,t2/2]);
nurbsArcTop = scaleNURBS(nurbsArcTop,t1/2);
nurbsArcTop = rotateNURBS(nurbsArcTop,'theta',pi/2,'rotAxis','Xaxis');
nurbsArcTop = translateNURBS(nurbsArcTop,[mean([R_8i,R_8o]),0,h9]);
nurbs = {nurbs{1}, nurbsArcBot{1}, nurbsArcTop{1}};

nurbs = revolveNURBS(nurbs);

nurbs = nurbs([2,1,3]);
nurbs([1,3]) = permuteNURBS(nurbs([1,3]),[3,2,1]);
nurbs(2) = permuteNURBS(nurbs(2),[3,1,2]);
if options.uniformDegree
    nurbs = makeUniformNURBSDegree(nurbs);
end





