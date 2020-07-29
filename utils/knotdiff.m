function d = knotdiff(Xi_1,Xi_2)

d = [];
uniqueXi_1 = unique(Xi_1);

for i = 1:length(uniqueXi_1)
    xi = uniqueXi_1(i);
    
    mult_1 = length(find(Xi_1 == xi));
    mult_2 = length(find(Xi_2 == xi));
    
    d = [d xi*ones(1,mult_1-mult_2)];
    
end