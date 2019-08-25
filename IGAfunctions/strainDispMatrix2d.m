function B = strainDispMatrix2d(nen,dRdx)

B(1,1:nen)          = dRdx(1,:);
B(2,nen+1:2*nen)    = dRdx(2,:);
B(3,1:nen)          = dRdx(1,:);
B(3,nen+1:2*nen)	= dRdx(2,:);