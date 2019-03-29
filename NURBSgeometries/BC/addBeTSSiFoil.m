function nurbs = addBeTSSiFoil(l_ls, b_ls, l_us, b_us, h_s, delta_s, x_0, alpha, exact)

if exact
    t = b_ls/l_ls;
    p = 8;
    b_u = getNACAapprox(t,p);
    t = b_us/l_us;
    [b_o, Xi] = getNACAapprox(t,p);
else
    N1 = 1;
    N2 = 7;
    p = 2;
    xi_T = 0.168967783470083;
    t = b_ls/l_ls;
    b_u = getNACAapprox(t,p,N1,N2,xi_T);
    t = b_us/l_us;
    [b_o, Xi] = getNACAapprox(t,p,N1,N2,xi_T);
end
controlPts = ones(4,p+1,4);
controlPts(1,:,1) = b_u(:,1)*l_ls;
controlPts(2,:,1) = b_u(:,2)*l_ls;
controlPts(3,:,1) = 0;
controlPts(1,:,2) = delta_s + b_o(:,1)*l_us;
controlPts(2,:,2) = b_o(:,2)*l_us;
controlPts(3,:,2) = h_s;
controlPts(1,:,:) = -controlPts(1,:,:);
controlPts(:,:,3) = controlPts(:,:,2);
controlPts(2,:,3) = -controlPts(2,:,3);
controlPts(:,:,4) = controlPts(:,:,1);
controlPts(2,:,4) = -controlPts(2,:,4);
controlPts = transformPts(controlPts,alpha,x_0);
Eta = [0,0,1,2,3,3]/3;

nurbs = createNURBSobject(controlPts,{Xi,Eta});
% nurbs = elevateNURBSdegree(nurbs,[0 1]);
function controlPts = transformPts(controlPts,alpha,x_0)

R_x = rotationXaxis(alpha);
ss = size(controlPts);
for i = 1:ss(2)
    for j = 1:ss(3)
        controlPts(1:3,i,j) = x_0.' + R_x*controlPts(1:3,i,j);
    end
end

function R_x = rotationXaxis(alpha)

R_x = [1, 0,           0;
       0, cos(alpha), -sin(alpha);
       0, sin(alpha),  cos(alpha)];
