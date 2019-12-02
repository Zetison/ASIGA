function nurbs = getShelveBracketData(r,R,L,l,h,t,s,d,h_front)

nurbs1 = getScrewData(r,R,l,h,t,2,s);
nurbs1 = rotateNURBS(nurbs1,-pi/2,'Zaxis');
nurbs1 = rotateNURBS(nurbs1,-pi/2,'Xaxis');
nurbs1 = translateNURBS(nurbs1,[0,0,l]);

nurbs2 = getScrewData(r,R,l,h,t,1,s);
nurbs2 = rotateNURBS(nurbs2,-pi/2,'Xaxis');
nurbs2 = translateNURBS(nurbs2,[0,0,l]);

nurbs(1) = nurbs1;
nurbs(2) = nurbs2;
nurbs(3:4) = translateNURBS(mirrorNURBS(nurbs([2,1]),'z'),[0,0,-L]);
controlPts = nurbs{3}.coeffs(:,3:end,end,:);
controlPts(:,:,2,:) = nurbs{2}.coeffs(:,3:end,end,:);
Xi = [0,0,0,1,1,1];
Eta = [0,0,1,1];
nurbs{5} = createNURBSobject(controlPts,{Xi,Eta, nurbs{3}.knots{3}});
nurbs{7} = flipNURBSparametrization(getQuarterCylinderData(l,h_front),'eta');
nurbs(7) = translateNURBS(nurbs(7),[0,t+d,-h_front]);


Xi = [0 0 0.5 1 1];
Eta = [0 0 1 1];
Zeta = [0 0 1 1];
controlPts = zeros(4,3,2,2);
q = L-h_front;
z0 = -q/2;  % free parameter
y0 = d/2;   % free parameter
% a = -1/(z0/y0); % free parameter
a = 0.5; % free parameter
b = z0-a*y0;
z1 = a*0 + b;
y2 = d/2; % free parameter
y1 = 4*d/5; % free parameter
z2 = a*y1 + b;

controlPts(:,1,1,1) = [ 0   d     0      1  	];
controlPts(:,2,1,1) = [ 0   y0    z0     1  	];
controlPts(:,3,1,1) = [ 0   0     -q     1  	];

controlPts(:,1,2,1) = [ 0   d     0      1  	];
controlPts(:,2,2,1) = [ 0   0     0      1  	];
controlPts(:,3,2,1) = [ 0   0     -q     1  	];

controlPts(:,1,1,2) = [ l   d     0      1  	];
controlPts(:,2,1,2) = [ l   y0    z0     1  	];
controlPts(:,3,1,2) = [ l   0     -q     1  	];

controlPts(:,1,2,2) = [ l   d     0      1  	];
controlPts(:,2,2,2) = [ l   0     0      1  	];
controlPts(:,3,2,2) = [ l   0     -q     1  	];

nurbs(8) = translateNURBS(createNURBSobject(controlPts,{Xi, Eta, Zeta}),[0,t,-h_front]);

nurbs([1:5,7]) = elevateDegreeInPatches(nurbs([1:5,7]),[0,1,1]);
nurbs(8) = elevateDegreeInPatches(nurbs(8),[1,1,1]);

controlPts = zeros(4,3,2,2);
Xi = [0 0 0 1 1 1];
Eta = [0 0 1 1];
Zeta = [0 0 1 1];

controlPts(:,1,1,1) = [ 0   d     -h_front      1  	];
controlPts(:,2,1,1) = [ 0   y2    -h_front      1  	];
controlPts(:,3,1,1) = [ 0   0     -h_front      1  	];

controlPts(:,1,2,1) = [ l   d     -h_front      1  	];
controlPts(:,2,2,1) = [ l   y2    -h_front      1  	];
controlPts(:,3,2,1) = [ l   0     -h_front      1  	];

controlPts(:,1,1,2) = [ 0   d     0      1  	];
controlPts(:,2,1,2) = [ 0   y2    0      1  	];
controlPts(:,3,1,2) = [ 0   0     0      1  	];

controlPts(:,1,2,2) = [ l   d     0      1  	];
controlPts(:,2,2,2) = [ l   y2    0      1  	];
controlPts(:,3,2,2) = [ l   0     0      1  	];

nurbs(6) = translateNURBS(createNURBSobject(controlPts,{Xi,Eta,Zeta}),[0,t,0]);

nurbs(6) = elevateDegreeInPatches(nurbs(6),[0,1,1]);
% nurbs(8:14) = mirrorNURBS(nurbs,'x');
