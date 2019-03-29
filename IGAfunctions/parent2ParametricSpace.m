function xi = parent2ParametricSpace(xiE,xiTilde)

xi = (xiTilde*(xiE(2) - xiE(1)) + (xiE(2) + xiE(1)))*0.5; % = xiE(1) + (xiTilde+1)*0.5*(xiE(2) - xiE(1)); 
