function B = strainDispMatrix(nen,dRdx)
d_p = size(dRdx,1);
switch d_p
    case 2
        B(1,1:nen)          = dRdx(1,:);
        B(2,nen+1:2*nen)    = dRdx(2,:);
        B(3,1:nen)          = dRdx(2,:);
        B(3,nen+1:2*nen)	= dRdx(1,:);
    case 3
        B(1,1:nen)          = dRdx(1,:);
        B(2,nen+1:2*nen)    = dRdx(2,:);
        B(3,2*nen+1:3*nen)  = dRdx(3,:);

        B(4,nen+1:2*nen)	= dRdx(3,:);
        B(4,2*nen+1:3*nen)  = dRdx(2,:);

        B(5,1:nen)          = dRdx(3,:);
        B(5,2*nen+1:3*nen)	= dRdx(1,:);

        B(6,1:nen)          = dRdx(2,:);
        B(6,nen+1:2*nen)    = dRdx(1,:);
end