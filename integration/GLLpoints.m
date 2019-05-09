function [points, weights] = GLLpoints(n)
points = zeros(1,n); 
weights = zeros(1,n);

if n > 64
	warning('There exist no listed gausspoints of order 65 or more: using 64 quadrature points instead...')
	n = 64;
end

switch n
	case 2
		% Set Gauss quadrature points
		points(2) =  1;
		
		% Set Gauss quadrature weights
		weights(2) = 1;
		
	case 3
		% Set Gauss quadrature points
