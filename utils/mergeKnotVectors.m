function A = mergeKnotVectors(Xi_1, Xi_2)

A = Xi_1;
for i = 1:length(Xi_2)
    xi = Xi_2(i);
    
    mult_1 = length(find(A == xi));
    mult_2 = length(find(Xi_2 == xi));
    
    A = [A xi*ones(1,mult_2 - mult_1)];
    A = sort(A);
        
end