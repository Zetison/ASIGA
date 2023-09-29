function nurbsCol = getBCinteriorData()


BCP = setBCParameters;
a = BCP.a;
b = BCP.b;
c = BCP.c;
s = BCP.s;

%% 6.1.8 Pressure hull
nurbsCol = cell(1,28);
layer = setPHParameters;
nurbsCol(1) = translateNURBS(subNURBS(getBeTSSiPHData,'at',[0,0;0,0;0,1]), [-(9.0-a)-layer{1}.gd/2,0,0]);

%% 6.1.9 Torpedo tubes
counter = 2;
for i = 0:3
    for j = [-1,1]
        r1 = 7.0;
        x0 = [BCP.a-5.7-r1/2, -1.5+i, j*0.5];
        nurbsCol(counter) = getCylinderData('L', r1, 'R',0.3, 'd_p',2, 'x_0', x0, 'alignWithAxis','Xaxis');
        counter = counter + 1;
    end
end
   
%% 6.1.10 Bulkheads
% Area 1
controlPts = zeros(4,2,2);
controlPts(:,1,1) = [9, 0, b*cos(30*pi/180), 1];
controlPts(:,2,1) = [9, 0, -1.2,             1];
controlPts(:,1,2) = [0, 0, b*cos(30*pi/180), 1];
controlPts(:,2,2) = [0, 0, -1.2,             1];
nurbsCol(10) = createNURBSobject(controlPts,{[0,0,1,1],[0,0,1,1]});


% Area 2
controlPts(:,1,1) = [5.5,  b, b*cos(30*pi/180), 1];
controlPts(:,2,1) = [5.5,  b, -1.2,             1];
controlPts(:,1,2) = [5.5, -b, b*cos(30*pi/180),1];
controlPts(:,2,2) = [5.5, -b, -1.2,            1];
nurbsCol(11) = createNURBSobject(controlPts,{[0,0,1,1],[0,0,1,1]});

% Area 3
controlPts(:,1,1) = [9,  b, b*cos(30*pi/180),  1];
controlPts(:,2,1) = [9, -b, b*cos(30*pi/180), 1];
controlPts(:,1,2) = [0,  b, b*cos(30*pi/180),  1];
controlPts(:,2,2) = [0, -b, b*cos(30*pi/180), 1];
nurbsCol(12) = createNURBSobject(controlPts,{[0,0,1,1],[0,0,1,1]});

% Area 4
controlPts(:,1,1) = [9,  b, -1.2,  1];
controlPts(:,2,1) = [9, -b, -1.2, 1];
controlPts(:,1,2) = [0,  b, -1.2,  1];
controlPts(:,2,2) = [0, -b, -1.2, 1];
nurbsCol(13) = createNURBSobject(controlPts,{[0,0,1,1],[0,0,1,1]});

% Area 5
controlPts(:,1,1) = [2.7, b, b*cos(30*pi/180), 1];
controlPts(:,2,1) = [2.7, b, -3.5,             1];
controlPts(:,1,2) = [2.7, -b, b*cos(30*pi/180),1];
controlPts(:,2,2) = [2.7, -b, -3.5,            1];
nurbsCol(14) = createNURBSobject(controlPts,{[0,0,1,1],[0,0,1,1]});
for i = 10:14
    nurbsCol(i) = translateNURBS(nurbsCol{i}, [-2,0,0]); % correct the location of the origin
end

%% 6.1.11 Interior structures of sail

% Area 6
controlPts(:,1,1) = [a-22, s, c,     1];
controlPts(:,2,1) = [a-22, -s, c,    1];
controlPts(:,1,2) = [a-22, s, c+3.5, 1];
controlPts(:,2,2) = [a-22, -s, c+3.5,1];
nurbsCol(15) = createNURBSobject(controlPts,{[0,0,1,1],[0,0,1,1]});

% Area 7
controlPts(:,1,1) = [a-25, s, c,     1];
controlPts(:,2,1) = [a-25, -s, c,    1];
controlPts(:,1,2) = [a-25, s, c+3.5, 1];
controlPts(:,2,2) = [a-25, -s, c+3.5,1];
nurbsCol(16) = createNURBSobject(controlPts,{[0,0,1,1],[0,0,1,1]});

% Area 8
controlPts(:,1,1) = [a-19, s, c+2.5, 1];
controlPts(:,2,1) = [a-19, -s,c+2.5, 1];
controlPts(:,1,2) = [a-25, s,c+2.5,  1];
controlPts(:,2,2) = [a-25, -s,c+2.5, 1];
nurbsCol(17) = createNURBSobject(controlPts,{[0,0,1,1],[0,0,1,1]});

% Left tube ø1.2m
ytop = 6.5;
ybot = 3;
z = 0;
x0 = [a-21, z, ybot];
nurbsCol(18) = getCylinderData('R', 1.2/2, 'L', ytop-ybot, 'd_p',2, 'x_0',x0, 'alignWithAxis','Zaxis');

% Middle tube ø0.6m
ytop = 7.5;
ybot = 3;
z = 0;
x0 = [a-22.8, z, ybot];
nurbsCol(19) = getCylinderData('R', 0.6/2, 'L', ytop-ybot, 'd_p',2, 'x_0',x0, 'alignWithAxis','Zaxis');

% Righ_s tubes ø0.2m
counter = 20;
ytop = 7.5;
ybot = 3;
for i = [-1,1]
    z = i*0.3;
    x0 = [a-23.2, z, ybot];
    nurbsCol(counter) = getCylinderData('R', 0.2/2, 'L', ytop-ybot, 'd_p',2, 'x_0',x0, 'alignWithAxis','Zaxis');
    counter = counter + 1;
end

%% 6.1.12 Interior of aft section
g = BCP.L;
% Area 9
controlPts(:,1,1) = [-g,     b*sin(30*pi/180), b*cos(30*pi/180), 1];
controlPts(:,2,1) = [-g,    -b*sin(30*pi/180), b*cos(30*pi/180), 1];
controlPts(:,1,2) = [-(g+3), b*sin(30*pi/180), b*cos(30*pi/180), 1];
controlPts(:,2,2) = [-(g+3),-b*sin(30*pi/180), b*cos(30*pi/180), 1];
nurbsCol(22) = createNURBSobject(controlPts,{[0,0,1,1],[0,0,1,1]});

% Area 10
controlPts(:,1,1) = [-g,        b*sin(30*pi/180), b*cos(30*pi/180), 1];
controlPts(:,2,1) = [-g,        b*sin(30*pi/180), tan(30*pi/180)*b*cos(30*pi/180), 1];
controlPts(:,1,2) = [-(g+3),    b*sin(30*pi/180), b*cos(30*pi/180), 1];
controlPts(:,2,2) = [-(g+3),    b*sin(30*pi/180), tan(30*pi/180)*b*cos(30*pi/180), 1];
nurbsCol(23) = createNURBSobject(controlPts,{[0,0,1,1],[0,0,1,1]});

% Area 11
controlPts(:,1,1) = [-g,        -b*sin(30*pi/180), b*cos(30*pi/180), 1];
controlPts(:,2,1) = [-g,        -b*sin(30*pi/180), tan(30*pi/180)*b*cos(30*pi/180), 1];
controlPts(:,1,2) = [-(g+3),    -b*sin(30*pi/180), b*cos(30*pi/180), 1];
controlPts(:,2,2) = [-(g+3),    -b*sin(30*pi/180), tan(30*pi/180)*b*cos(30*pi/180), 1];
nurbsCol(24) = createNURBSobject(controlPts,{[0,0,1,1],[0,0,1,1]});

% Area 12
controlPts(:,1,1) = [-(g+3), -b*sin(30*pi/180), b*cos(30*pi/180), 1];
controlPts(:,2,1) = [-(g+3), -b*sin(30*pi/180), tan(30*pi/180)*b*cos(30*pi/180), 1];
controlPts(:,1,2) = [-(g+3),  b*sin(30*pi/180), b*cos(30*pi/180), 1];
controlPts(:,2,2) = [-(g+3),  b*sin(30*pi/180), tan(30*pi/180)*b*cos(30*pi/180), 1];
nurbsCol(25) = createNURBSobject(controlPts,{[0,0,1,1],[0,0,1,1]});

%% 6.1.12 Sonar-Array and Sonar-Window
% Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;

% nurbs = rotateNURBS(getDiskData('R',3.6/2, 't', 0, 'Xi',Xi),'theta',pi/2,'rotAxis','Yaxis');
% counter = 26;
% for i = [1,-1]
%     x0 = [a-3.2, 0, -1.7+i*0.6/2];
%     nurbsCol(counter) = translateNURBS(nurbs{1}, x0);
%     counter = counter + 1;
% end
x0 = [a-3.2, 0, -1.7-0.6/2];
nurbsCol(26:28) = subNURBS(getCylinderData('R', [3.6/2, 0], 'L', 0.6, 'd_p',3, 'x_0', x0, 'alignWithAxis','Zaxis'),'at',[0,1;0,0;1,1]);

