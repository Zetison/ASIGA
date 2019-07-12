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
controlPts(:,1,1) = [ R/2   0       0      1            ];
controlPts(:,2,1) = [ R/2   R/4 	0      w            ];
controlPts(:,3,1) = [ R/2   R/2     0      1            ];
controlPts(:,1,2) = [ R/4   0       0      w            ];
controlPts(:,2,2) = [ R/4   R/4 	0      w^2          ];
controlPts(:,3,2) = [ R/4   R/2     0      w            ];
controlPts(:,1,3) = [ 0     0       0      1            ];
controlPts(:,2,3) = [ 0     R/4 	0      w            ];
controlPts(:,3,3) = [ 0     R/2     0      1            ];

nurbs{1} = createNURBSobject(controlPts,{Xi, Eta});

v1 = [R, (sqrt(2)-1)*R];
v2 = [R/2,R/4];
v3 = (v1+v2)/2;
v4 = [R/sqrt(2),R/sqrt(2)];
v5 = (v4+[R/2,R/2])/2;
controlPts(:,1,1) = [ R             0       0      1            ];
controlPts(:,2,1) = [ v1                    0      w            ];
controlPts(:,3,1) = [ v4                    0      1            ];
controlPts(:,1,2) = [ 3*R/4         0       0      1            ];
controlPts(:,2,2) = [ v3                    0      w            ];
controlPts(:,3,2) = [ v5                    0      1            ];
controlPts(:,1,3) = [ R/2           0       0      1            ];
controlPts(:,2,3) = [ R/2           R/4 	0      w            ];
controlPts(:,3,3) = [ R/2           R/2     0      1            ];

nurbs{2} = createNURBSobject(controlPts,{Xi, Eta});

v1 = [(sqrt(2)-1)*R, R];
v2 = [R/4, R/2];
v3 = (v1+v2)/2;
v5 = ([R/sqrt(2),R/sqrt(2)]+[R/2,R/2])/2;
controlPts(:,1,1) = [ v4                0      1            ];
controlPts(:,2,1) = [ v1                0      w            ];
controlPts(:,3,1) = [ 0  	R           0      1            ];
controlPts(:,1,2) = [ v5                0      1            ];
controlPts(:,2,2) = [ v3                0      w            ];
controlPts(:,3,2) = [ 0 	3*R/4       0      1            ];
controlPts(:,1,3) = [ R/2 	R/2         0      1            ];
controlPts(:,2,3) = [ R/4	R/2         0      w            ];
controlPts(:,3,3) = [ 0 	R/2         0      1            ];

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

