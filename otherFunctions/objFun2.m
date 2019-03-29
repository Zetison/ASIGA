function fun = objFun2(nurbs, xi,eta)
X = evaluateNURBS(nurbs,[xi, eta]);
fun = abs(X(3));