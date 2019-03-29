function controlPts = calcTensorRotCtrlPts(ctrlPtsXi,ctrlPtsEta)
n_xi = size(ctrlPtsXi,2);
n_eta = size(ctrlPtsEta,2);
controlPts = zeros(4,n_xi,n_eta);
for i = 1:n_xi
    for j = 1:n_eta
        controlPts(4,i,j) = ctrlPtsXi(3,i)*ctrlPtsEta(3,j);
        controlPts(1,i,j) = ctrlPtsEta(1,j);
        controlPts(2:3,i,j) = ctrlPtsXi(1:2,i)*ctrlPtsEta(2,j);
    end
end
