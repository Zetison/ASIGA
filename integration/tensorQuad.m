function [Q, W] = tensorQuad(uOrder,vOrder,wOrder)

switch nargin
    case 1
        [Q, W] = getQuadFromFile(uOrder);

    case 2
        Q  = zeros(uOrder*vOrder,2);
        W = zeros(uOrder*vOrder,1);

        % Gauss quadrature along two directions

        [ptU, wtU] = getQuadFromFile(uOrder);
        [ptV, wtV] = getQuadFromFile(vOrder);

        A = 1;

        for j = 1:vOrder
            for i = 1:uOrder
                Q(A,:) = [ptU(i), ptV(j)];
                W(A)  = wtU(i)*wtV(j);
                A = A+1;
            end
        end
    case 3
        Q = zeros(uOrder*vOrder*wOrder,3);
        W = zeros(uOrder*vOrder*wOrder,1);

        % Gauss quadrature along three directions

        [ptU, wtU] = getQuadFromFile(uOrder);
        [ptV, wtV] = getQuadFromFile(vOrder);
        [ptW, wtW] = getQuadFromFile(wOrder);

        A = 1;

        for k = 1:wOrder
            for j = 1:vOrder
                for i = 1:uOrder
                    Q(A,:) = [ptU(i), ptV(j), ptW(k)];
                    W(A)  = wtU(i)*wtV(j)*wtW(k);
                    A = A+1;
                end
            end
        end
end


