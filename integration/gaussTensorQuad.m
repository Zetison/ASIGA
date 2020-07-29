function [Q, W] = gaussTensorQuad(order,type)
if nargin < 2
    type = 'Legendre';
end
switch numel(order)
    case 1
        [Q, W] = eval(['gauss' type 'Quad(order(1))']);

    case 2
        quadpoint  = zeros(order(1)*order(2),2);
        quadweight = zeros(order(1)*order(2),1);

        % Gauss quadrature along two directions

        [ptU, wtU] = eval(['gauss' type 'Quad(order(1))']);
        [ptV, wtV] = eval(['gauss' type 'Quad(order(2))']);

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
    case 3
        quadpoint  = zeros(order(1)*order(2)*order(3),3);
        quadweight = zeros(order(1)*order(2)*order(3),1);

        % Gauss quadrature along three directions

        [ptU, wtU] = eval(['gauss' type 'Quad(order(1))']);
        [ptV, wtV] = eval(['gauss' type 'Quad(order(2))']);
        [ptW, wtW] = eval(['gauss' type 'Quad(order(3))']);

        A = 1;

        for k = 1:order(3)
            for j = 1:order(2)
                for i = 1:order(1)
                    quadpoint(A,:) = [ptU(i), ptV(j), ptW(k)];
                    quadweight(A)  = wtU(i)*wtV(j)*wtW(k);
                    A = A+1;
                end
            end
        end

        Q = quadpoint;
        W = quadweight;
end

