function [nodes, visElements] = createMesh(varCol, xBnd, zBnd, Nx, Nz)

Eta = varCol.knots{2};

Xi = [Eta (1-fliplr(Eta))];

nodes = zeros(Nx*Nz, 3);
visElements = zeros((Nx-1)*(Ny-1),4);
x_arr = linspace(xBnd(1),xBnd(2),Nx);
z_arr = linspace(zBnd(1),zBnd(2),Nz);
counter = 1;
eVis = 1;
nurbs = varCol.nurbs;
for j = 1:Nz
    z = z_arr(j);
    for i = 1:Nx
        x = x_arr(i);
        X = [x, 0, z];
        
        xi = 0;
        objfun = @(eta) objfunction(xi, eta);
        eta = bisection(objfun, 0, 1, 1000, 1e-10);
        P1 = evaluateNURBS(nurbs, [xi, eta, 0]);
        
        xi = 0.5;
        objfun = @(eta) objfunction(xi, eta);
        eta = bisection(objfun, 0, 1, 1000, 1e-10);
        P2 = evaluateNURBS(nurbs, [xi, eta, 0]);
        if norm(X-P1) < norm(X-P2)
            P = P1;
        else
            P = P2;
        end
        [~, ~, dCdeta, ~] = evaluateNURBS_deriv(nurbs, [xi, eta, 0]);
        n = dCdeta/norm(dCdeta);
        if dot(n,
        
        
        if r == 1
            nodes(counter) = X;
            counter = counter + 1;
            if i ~= 1 && j ~= 1
                visElements(eVis,1) = i   +   (j-1)*noXiKnots;
                visElements(eVis,2) = i+1 +   (j-1)*noXiKnots;
                visElements(eVis,3) = i+1 +       j*noXiKnots;
                visElements(eVis,4) = i   +       j*noXiKnots;
                eVis = eVis + 1;
            end
        end
    end
end


function f = objfunction(X, xi, eta)

[C, ~, dCdeta, ~] = evaluateNURBS_deriv(nurbs, [xi, eta, 0]);
f = dot(X-C, dCdeta);



