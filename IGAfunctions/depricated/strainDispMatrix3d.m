function B = strainDispMatrix3d(nen,dRdx)
error('Depricated, use strainDispMatrix(nen,dRdx) instead')

B(1,1:nen)          = dRdx(1,:);
B(2,nen+1:2*nen)    = dRdx(2,:);
B(3,2*nen+1:3*nen)  = dRdx(3,:);

B(4,nen+1:2*nen)	= dRdx(3,:);
B(4,2*nen+1:3*nen)  = dRdx(2,:);

B(5,1:nen)          = dRdx(3,:);
B(5,2*nen+1:3*nen)	= dRdx(1,:);

B(6,1:nen)          = dRdx(2,:);
B(6,nen+1:2*nen)    = dRdx(1,:);