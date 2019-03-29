function nurbs = insertKnotsInNURBS_NonLinearParametrization(nurbs, noNewPoints)


if strcmp(nurbs.type, '3Dvolume')
    
    
elseif strcmp(nurbs.type, '2Dsurface')
    Xi = nurbs.knots{1};
    uniqueXi = unique(Xi);
    for i = 1:length(uniqueXi)-1
        
        xi_1 = uniqueXi(i);
        xi_2 = uniqueXi(i+1);
        P1_xi_1 = evaluateNURBS(nurbs, [xi_1, 0]);
        P1_xi_2 = evaluateNURBS(nurbs, [xi_2, 0]);
        
        P1 = (P1_xi_1 + P1_xi_2)/2;
        
        P2_xi_1 = evaluateNURBS(nurbs, [xi_1, 1]);
        P2_xi_2 = evaluateNURBS(nurbs, [xi_2, 1]);
        
        P2 = (P2_xi_1 + P2_xi_2)/2;
        
        dist1 = @(xi) dot(P1 - evaluateNURBS(nurbs, [xi, 0]), evaluateNURBS_deriv2(nurbs, [xi 0], 'xi'));
        
        xi_1_bar = bisection(dist1, xi_1, xi_2, 100, 1e-10);
        
        dist2 = @(xi) dot(P2 - evaluateNURBS(nurbs, [xi, 1]), evaluateNURBS_deriv2(nurbs, [xi 1], 'xi'));
        
        xi_2_bar = bisection(dist2, xi_1, xi_2, 100, 1e-10);
        
        norm1 = norm(evaluateNURBS_deriv2(nurbs, [xi_1_bar, 0], 'xi'));
        norm2 = norm(evaluateNURBS_deriv2(nurbs, [xi_2_bar, 1], 'xi'));

        xi_bar = (norm1*xi_1_bar + norm2*xi_2_bar)/(norm1 + norm2);
        
        nurbs = insertKnotsInNURBS(nurbs,{xi_bar []});
    end   
elseif strcmp(nurbs.type, '2Dcurve')
    
    
elseif strcmp(nurbs.type, '1Dnurbs')
    for i = 1:length(Xi)-1
        xi_1 = Xi(i);
        xi_2 = Xi(i+1);
        P_xi_1 = evaluateNURBS(nurbs, xi_1);
        P_xi_2 = evaluateNURBS(nurbs, xi_2);
         
        Xi_new = zeros(noNewPoints,1);
        for j = 1:noNewPoints
            P_j = P_xi_1 + j/(noNewPoints+1)*(P_xi_1 + P_xi_2);
            Xi_new(j) = pointInversion(nurbs,P_j,1e-10);
        end
        nurbs = insertKnotsInNURBS(nurbs,Xi_new);
    end
end