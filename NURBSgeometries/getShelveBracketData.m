function nurbs = getShelveBracketData(r,R,L,l,h,t,s,d)

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
controlPts = nurbs{5}.coeffs(:,:,:,end);
controlPts(:,:,:,2) = controlPts;
controlPts(2,:,:,2) = controlPts(2,:,:,2) + d;
controlPts = permute(controlPts,[1,4,2,3]); % done in order to get zeta = 1 on top
nurbs{6} = createNURBSobject(controlPts,{Eta,Xi,Eta});
nurbs(1:5) = elevateDegreeInPatches(nurbs(1:5),[0,1,1]);
nurbs(6) = elevateDegreeInPatches(nurbs(6),[1,0,1]);
% nurbs(7:12) = mirrorNURBS(nurbs,'x');
