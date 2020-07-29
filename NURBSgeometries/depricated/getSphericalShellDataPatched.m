function nurbs = getSphericalShellDataPatched(R_o,t)
error('Use getEllipsoidData() instead')

nurbs = cell(1,6);
nurbs(1) = getSphericalShellPatchData(R_o,t);

for i = 2:4
    nurbs(i) = rotateNURBS(nurbs{1},(i-1)*pi/2,'Yaxis');
end
nurbs(5) = rotateNURBS(nurbs{1},pi/2,'Xaxis');
nurbs(6) = rotateNURBS(nurbs{1},-pi/2,'Xaxis');
