function nurbs = getModel3Data(R_1, R_2, t, L, nElXi)
if nargin < 5
    nElXi = 4;
end
n_xi = 2*nElXi+1;
Xi = zeros(1,n_xi+3);
Xi(end-2:end) = 1;
for i = 1:nElXi
    Xi(2*i+2:2*i+3) = i/nElXi;
end
totLength = R_1*pi/2 + sqrt(L^2+(R_1-R_2)^2)+ R_2*pi/2;
Eta = zeros(1, 10);
Eta(4:5) = R_1*pi/2/totLength;
Eta(6:7) =(R_1*pi/2 + sqrt(L^2+(R_1-R_2)^2))/totLength;
Eta(8:end) = 1;
Zeta = [0 0 1 1];

controlPts = zeros(4,n_xi,7);
controlPts(:,:,:,1) = getCtrlPts(R_1-t, R_2-t, L, Xi);
controlPts(:,:,:,2) = getCtrlPts(R_1, R_2, L, Xi);
nurbs = createNURBSobject(controlPts,{Xi, Eta, Zeta});


function controlPts = getCtrlPts(R_1, R_2, L, Xi)

R_m = (R_1+R_2)/2;

ctrlPtsXi = parmArc(Xi,2*pi);
ctrlPtsEta = [-L-R_1   0    1;
              -L-R_1   R_1	1/sqrt(2);
              -L       R_1  1;
              -L/2     R_m  1;
               0       R_2  1;
               R_2     R_2 	1/sqrt(2);
               R_2     0    1].';

controlPts = calcTensorRotCtrlPts(ctrlPtsXi,ctrlPtsEta);

