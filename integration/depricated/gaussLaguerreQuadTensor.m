function [W,Q] = gaussLaguerreQuadTensor(order, r)
error('use gaussTensorQuad() instead')
switch nargin
    case 2
        quadpoint  = zeros(order(1)*order(2),2);
        quadweight = zeros(order(1)*order(2),1);

        % Gauss quadrature along two directions

        [ptU, wtU] = gaussLaguerreQuad(order(1),r(1));
        [ptV, wtV] = gaussLaguerreQuad(order(2),r(2));

        A = 1;

        for j = 1:order(2)
            for i = 1:order(1)
                quadpoint(A,:) = [ptU(i), ptV(j)];
                quadweight(A)  = wtU(i)*wtV(j);
                A = A+1;
            end
        end

        Q = quadpoint;
        W = quadweight;
end

