function nurbs = addSailToBeTSSi(nurbs)


eta1 = 2/3;
eta2 = 5/6;
lift = 2;

nurbs = insertKnotsInNURBS(nurbs,{0.5 linspace2(eta1,eta2,7)});

idxXi = 14:15;
idxEta = 12:13;
nurbs.coeffs(3,idxXi,idxEta) = nurbs.coeffs(3,idxXi,idxEta) + lift;
nurbs.coeffs(1,idxXi,idxEta(2)) = nurbs.coeffs(1,idxXi,idxEta(2)+1);
nurbs.coeffs(2,idxXi,idxEta(1)) = 0.2*nurbs.coeffs(2,idxXi,idxEta(1));
nurbs.coeffs(2,idxXi,idxEta(2)) = 2*nurbs.coeffs(2,idxXi,idxEta(2));