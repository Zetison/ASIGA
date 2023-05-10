function nurbs = addBeTSSiFoil3(l_ls, b_ls, l_us, b_us,x_s,z,h_s,dx,dy,Xi,idx,weigh_ss,s2, r_o, r_u, axisFlipIdx, zFlip)
axisFlipIdx2 = 2;
if ~iscell(r_o)
    r_o = {r_o,r_o};
end
if nargin < 16
    axisFlipIdx = 2;
    axisFlipIdx2 = 3;
    zFlip = 0;
end
p_xi = numel(find(Xi == 0))-1;
[s_rep, I] = getRepeatedKnots(Xi,p_xi);
s1 = s_rep(idx);
controlPts = [0, 1; 1, 1];
nurbs1D_1 = createNURBSobject(controlPts,[0,0,1,1]);
nurbs1D_1 = elevateNURBSdegree(nurbs1D_1,p_xi-1);
nurbs1D_1 = insertKnotsInNURBS(nurbs1D_1,{{Xi(and(Xi < s1-10*eps, Xi > 0))/s1}});
XiTemp = [zeros(1,p_xi+1), (Xi(and(Xi > (s1+10*eps), Xi < (s2+10*eps)))-s1)/(s2-s1), 1];
xiT1 = aveknt(XiTemp, p_xi+1);
X1 = @(xi) x_s + dx - xi*(s2*l_ls+dx);
XiTemp2 = [zeros(1,p_xi+1), (Xi(Xi > s2+10*eps)-s2)/(1-s2)];
xiT2 = aveknt(XiTemp2, p_xi+1);
X2 = @(xi) x_s - s2*l_ls - xi*(1-s2)*l_ls;

n = numel(Xi)-(p_xi+1);
if any(isnan(weigh_ss))
    weigh_ss = ones(1,n);
end
controlPts1 = ones(4,n,3,2);
for i = 1:2
    temp1 = X1(xiT1);
    controlPts1(1,I(idx):I(idx+1),1,i) = temp1;
    temp2 = X2(xiT2);
    controlPts1(1,I(idx+1):end,1,i) = temp2;
    if dy(i) == 0
        t = b_ls/l_ls;
        xValues = controlPts1(1,I(idx+2)-1:I(idx+2)+1,1,i).';
        b_u = getNACAapprox3(t,p_xi,Xi,r_u,weigh_ss,idx+2,xValues,[I(idx+2)-1, I(idx+2)+1], axisFlipIdx, axisFlipIdx2);
    else
        t = b_ls/l_ls;
        b_u = getNACAapprox3(t,p_xi,Xi,r_u,weigh_ss);
    end
    t = b_us/l_us;
    b_o = getNACAapprox3(t,p_xi,Xi,r_o{i},weigh_ss);
    if b_us == 2 % sail
        b_u(:,3) = 4;
        b_o(:,3) = 7.5;
    elseif b_us == 0.3 % main rudders
        b_o(:,3) = 3.5;
    else % depth rudders
        b_o(:,2) = 3.5;
    end
    
    y = dy(i)+b_ls/2;
    controlPts1(1,1:I(idx),1,i) = dx+x_s;
    controlPts1(axisFlipIdx,1:I(idx),1,i) = zFlip + nurbs1D_1{1}.coeffs(1,:)*y;
    controlPts1(axisFlipIdx,I(idx):end,1,i) = zFlip + y;
    controlPts1(axisFlipIdx2,:,1,i) = z;
    controlPts1(4,:,1,i) = weigh_ss;
    
    controlPts1(1:3,:,2,i) = b_u.';
    controlPts1(4,:,2,i) = weigh_ss;
    controlPts1(1:3,:,3,i) = b_o.';
    controlPts1(axisFlipIdx2,:,3,i) = z+h_s;
    controlPts1(4,:,3,i) = weigh_ss;
end
controlPts = ones(4,n,6);
controlPts(:,:,1:3) = controlPts1(:,:,:,1);
controlPts(:,:,4:6) = controlPts1(:,:,3:-1:1,2);
controlPts(axisFlipIdx,:,4:6) = 2*zFlip-controlPts(axisFlipIdx,:,4:6);

Eta = [0,0,1,2,3,4,5,5]/5;

nurbs = createNURBSobject(controlPts,{Xi,Eta});
