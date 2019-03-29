function nurbs = getRodData(leftStartingPoint,L)

Xi = [0 0 1 1];

controlPts = zeros(2,2);

controlPts(:,1) = [leftStartingPoint; 	 1];
controlPts(:,2) = [leftStartingPoint+L; 1];


nurbs = createNURBSobject(controlPts,Xi);