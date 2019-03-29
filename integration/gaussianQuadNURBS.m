function [W,Q] = gaussianQuadNURBS(uOrder,vOrder,wOrder)

switch nargin
    case 1
        [Q, W] = gaussQuad(uOrder);

    case 2
        quadpoint  = zeros(uOrder*vOrder,2);
        quadweight = zeros(uOrder*vOrder,1);

        % Gauss quadrature along two directions

        [ptU, wtU] = gaussQuad(uOrder);
        [ptV, wtV] = gaussQuad(vOrder);

        A = 1;

        for j = 1:vOrder
            for i = 1:uOrder
                quadpoint(A,:) = [ptU(i), ptV(j)];
                quadweight(A)  = wtU(i)*wtV(j);
                A = A+1;
            end
        end

        Q = quadpoint;
        W = quadweight;
    case 3
        quadpoint  = zeros(uOrder*vOrder*wOrder,3);
        quadweight = zeros(uOrder*vOrder*wOrder,1);

        % Gauss quadrature along three directions

        [ptU, wtU] = gaussQuad(uOrder);
        [ptV, wtV] = gaussQuad(vOrder);
        [ptW, wtW] = gaussQuad(wOrder);

        A = 1;

        for k = 1:wOrder
            for j = 1:vOrder
                for i = 1:uOrder
                    quadpoint(A,:) = [ptU(i), ptV(j), ptW(k)];
                    quadweight(A)  = wtU(i)*wtV(j)*wtW(k);
                    A = A+1;
                end
            end
        end

        Q = quadpoint;
        W = quadweight;
end

