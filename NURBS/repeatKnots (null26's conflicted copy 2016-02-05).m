function nurbs = repeatKnots(nurbs)


p = nurbs.degree(1);
q = nurbs.degree(2);
r = nurbs.degree(3);

Xi = nurbs.knots{1};
Eta = nurbs.knots{2};
Zeta = nurbs.knots{3};

for i = 1:length(Xi)
    mm = length(find(Xi == Xi(i)));    
    nurbs = insertKnotsInNURBS(nurbs,{Xi(i)*ones(p-mm,1) [] []});
end

for j = 1:length(Eta)
    mm = length(find(Eta == Eta(j)));    
    nurbs = insertKnotsInNURBS(nurbs,{[] Eta(j)*ones(q-mm,1) []});
end

for k = 1:length(Zeta)
    mm = length(find(Zeta == Zeta(k)));    
    nurbs = insertKnotsInNURBS(nurbs,{[] [] Zeta(k)*ones(r-mm,1)});
end