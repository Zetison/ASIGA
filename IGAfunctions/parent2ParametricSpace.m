function xi = parent2ParametricSpace(xiE,xiTilde)

if size(xiE,1) > 1
    xi = zeros(size(xiTilde));
    for i = 1:size(xiE,1)
        xi(:,i) = parent2ParametricSpace(xiE(i,:),  xiTilde(:,i));
    end
else
    xi = (xiTilde*(xiE(2) - xiE(1)) + (xiE(2) + xiE(1)))*0.5; % = xiE(1) + (xiTilde+1)*0.5*(xiE(2) - xiE(1)); 
end
