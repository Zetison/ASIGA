function nurbs = getSphericalShellDataPatchedQuarter(R)
% data based on data from Cobb1988tts
Xi = [0 0 0 0 0 1 1 1 1 1];
Eta = [0 0 0 0 0 1 1 1 1 1];

if isa(R,'sym')
    sr2 = sqrt(sym('2'));                                                        
    sr3 = sqrt(sym('3'));                                                        
    sr6 = sqrt(sym('6'));   
else
    sr2 = sqrt(2);                                                        
    sr3 = sqrt(3);                                                        
    sr6 = sqrt(6);   
end
controlPts = zeros(4,5,5,class(R)); 
w1 = (6 + 2*sr2 + sr3 + 2*sr6)/4;
w2 = (46 + 6*sr2 + sr3 + 16*sr6)/24;
w3 = (22 - 2*sr2 - 3*sr3 + 8*sr6)/8;
w4 = (-2*sr2 + 3*sr6 - 2*sr3 + 8)/2;
w5 = (166 - 11*sr3 + 60*sr6)/72;
w6 = (74 - 12*sr2 - 13*sr3 + 28*sr6)/24;
w7 = (-6*sr2 + 9*sr6 - 8*sr3 + 28)/6;
w8 = (-8*sr2 + 12*sr6 - 7*sr3 + 30)/8;
w9 = (-2*sr2 + 3*sr6 - 4*sr3 + 12)/2;
w10 = 4*(3-sr3);

a2 = (-2 + 4*sr2 + 3*sr3 - sr6)/8;
a3 = 2*a2;
a4 = (26 + 18*sr2 + 11*sr3 + 8*sr6)/24;
a5 = -3*(2 - 4*sr2 - 3*sr3 + sr6)/8;
a6 = (2 + 10*sr2 + 7*sr3)/8;
a7 = (3*sr3 + 4*sr2 - sr6 - 2)/2;

a8 = (-14 + 22*sr2 + 19*sr3 - 5*sr6)/48;
a9 = (-6 + 6*sr2 + 7*sr3 - sr6)/16;
a10 = (sr2 + 2*sr3 - 2)/4;
a11 = 2*a8;
a12 = (-22 + 34*sr2 + 29*sr3 - 8*sr6)/24;
a13 = (18 + 24*sr2 + 15*sr3 + 4*sr6)/24;
a14 = 2*a9;
a15 = (2 + 36*sr2 + 23*sr3 - 4*sr6)/24;
a16 = (10*sr3 + 12*sr2 - 3*sr6 - 8)/6;
a17 = 2*a10;
a18 = (-10 + 10*sr2 + 11*sr3 - 2*sr6)/8;
a19 = (9*sr3 + 16*sr2 - 4*sr6 - 2)/8;
a20 = (4*sr3 + 4*sr2 - sr6 - 4)/2;
a21 = (sr2 + 4*sr3 - 4)/2;
a22 = 4*(sr3-1);

controlPts(:,1,1) = [0  0  w1  w1];
controlPts(:,2,1) = [a2  0  w1  w1];
controlPts(:,3,1) = [a3  0  a4  w2];
controlPts(:,4,1) = [a5  0  a6  w3];
controlPts(:,5,1) = [a7  0  a7  w4];

controlPts(:,2,2) = [a2  a2  w1  w1];
controlPts(:,3,2) = [a3  a8  a4  w2];
controlPts(:,4,2) = [a5  a9  a6  w3];
controlPts(:,5,2) = [a7  a10  a7  w4];

controlPts(:,3,3) = [a11  a11  a13  w5];
controlPts(:,4,3) = [a12  a14  a15  w6];
controlPts(:,5,3) = [a16  a17  a16  w7];

controlPts(:,4,4) = [a18  a18  a19  w8];
controlPts(:,5,4) = [a20  a21  a20  w9];

controlPts(:,5,5) = [a22  a22  a22  w10];

for i = 1:5
    for j = i+1:5
        controlPts(1:2,i,j) = controlPts([2,1],j,i);
        controlPts(3:4,i,j) = controlPts(3:4,j,i);
    end
end  

for i = 1:5
    for j = 1:5
        controlPts(1:3,i,j) = controlPts(1:3,i,j)/controlPts(4,i,j);
    end
end  
controlPts(1:3,:,:) = controlPts(1:3,:,:)*R; 
nurbs = createNURBSobject(controlPts,{Xi, Eta}); 
