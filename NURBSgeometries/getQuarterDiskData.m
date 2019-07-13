function nurbs = getQuarterDiskData(R,alignWithAxis, x_0)
if nargin < 2
    alignWithAxis = 'Xaxis';
end
if nargin < 3
    x_0 = [0, 0, 0];
end

% create square
controlPts = zeros(4,3,3);
theta = pi/4;
w = cos(theta/2);
Xi = [0,0,0,1,1,1];
Eta = [0,0,0,1,1,1];
u1 = [0,0];
u2 = [R/2,0];
u3 = [R,0];
u4 = [0,R/2];
u5 = [1,1]*(sqrt(3)+3)/6*R/2;
u6 = [1,1]*R/sqrt(2);
u7 = [0,R];
v1 = [1,1]*(sqrt(3)+9)/12*R/4; % = ((u1+u4)/2+(u2+u5)/2)/2

controlPts(:,1,1) = [ u2            0      1            ];
controlPts(:,2,1) = [ (u2+u5)/2 	0      w            ];
controlPts(:,3,1) = [ u5            0      1            ];
controlPts(:,1,2) = [ (u1+u2)/2   	0      w            ];
controlPts(:,2,2) = [ v1            0      w^2          ];
controlPts(:,3,2) = [ (u4+u5)/2   	0      w            ];
controlPts(:,1,3) = [ u1            0      1            ];
controlPts(:,2,3) = [ (u1+u4)/2     0      w            ];
controlPts(:,3,3) = [ u4            0      1            ];

nurbs{1} = createNURBSobject(controlPts,{Xi, Eta});

v1 = R*[1, sqrt(2)-1];
controlPts(:,1,1) = [ u3                    0      1            ];
controlPts(:,2,1) = [ v1                    0      w            ];
controlPts(:,3,1) = [ u6                    0      1            ];
controlPts(:,1,2) = [ (u2+u3)/2             0      1            ];
controlPts(:,2,2) = [ ((u2+u5)/2+v1)/2      0      w            ];
controlPts(:,3,2) = [ (u5+u6)/2             0      1            ];
controlPts(:,1,3) = [ u2                    0      1            ];
controlPts(:,2,3) = [ (u2+u5)/2             0      w            ];
controlPts(:,3,3) = [ u5                    0      1            ];

nurbs{2} = createNURBSobject(controlPts,{Xi, Eta});

v1 = fliplr(v1);
controlPts(:,1,1) = [ u6                0      1            ];
controlPts(:,2,1) = [ v1                0      w            ];
controlPts(:,3,1) = [ u7                0      1            ];
controlPts(:,1,2) = [ (u5+u6)/2         0      1            ];
controlPts(:,2,2) = [ ((u4+u5)/2+v1)/2  0      w            ];
controlPts(:,3,2) = [ (u4+u7)/2         0      1            ];
controlPts(:,1,3) = [ u5                0      1            ];
controlPts(:,2,3) = [ (u4+u5)/2         0      w            ];
controlPts(:,3,3) = [ u4                0      1            ];

nurbs{3} = createNURBSobject(controlPts,{Xi, Eta});
switch alignWithAxis
    case 'Zaxis'
        % do nothing
    case 'Yaxis'
        nurbs = rotateNURBS(nurbs,-pi/2,'Zaxis');
        nurbs = rotateNURBS(nurbs,-pi/2,'Xaxis');
    case 'Xaxis'
        nurbs = rotateNURBS(nurbs,pi/2,'Zaxis');
        nurbs = rotateNURBS(nurbs,pi/2,'Yaxis');
end

nurbs = translateNURBS(nurbs,x_0);

