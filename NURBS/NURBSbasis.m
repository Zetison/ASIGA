function vargout = NURBSbasis(I,xi,degree,knots,weights,n)
% This routine compute the nonzero NURBS functions 
% and corresponding derivatives at pt
% This routine can be generalized based on Algorithm A4.2 in The NURBS Book Piegl1997tnb available at
% https://doi.org/10.1016/S0167-8396(98)00028-4

% Input
%       xi:                 evaluation points
%       degree:             NURBS degrees
%       knots:              knot vectors
%       weights:            NURBS weights

% Output
%       R, dR{1}, dR{2}, dR{3}:      array of the (p+1)(q+1)(r+1) NURBS 
%                                       functions and its derivatives which 
%                                       are nonzero at (xi,eta,zeta)
if nargin < 6
    n = 1;
end
d_p = numel(degree);
prec = class(xi);
noxi = numel(xi(:,1));
nen = prod(degree+1);
R = zeros(noxi,nen,1,prec);
W = zeros(noxi,1,1,prec);
B = cell(1,d_p);
if n > 0
    dW = cell(1,d_p);
    dR = cell(1,d_p);
else
    dR = {};
end
for i = 1:d_p
    B{i} = BsplineBasis(I(i), xi(:,i), degree(i), knots{i}, n);
    if n > 0
        dW{i} = zeros(noxi,1,n,prec);
        dR{i} = zeros(noxi,nen,n,prec);
    end
end
dd = ones(d_p);
for i = 1:d_p
    dd(i,i) = 2;
end
counter = 1;
switch d_p
    case 1
        W = B{1}(:,:,1)*weights;
        for i = 1:d_p
            for j = 1:n
                dW{i}(:,:,j) = B{i}(:,:,j+1,i)*weights;
            end
        end

        counter = 1;
        for k1 = 1:degree(1)+1 
            fact = weights(counter)./(W.*W);

            BB = B{1}(:,k1,1);
            R(:,counter) = BB.*fact.*W;

            if n > 0
                dWdxi = dW{i}(:,:,1);
                dR{1}(:,counter,1) = (B{1}(:,k1,2).*W - BB.*dWdxi).*fact;
            end
            if n > 1
                d2Wdxi2 = dW{i}(:,:,2);
                dR{1}(:,counter,2) = (B{1}(:,k1,3).*W - B{1}(:,k1,1).*d2Wdxi2).*fact - 2*fact./W.*(B{1}(:,k1,2).*W - B{1}(:,k1,1).*dWdxi).*dWdxi;
            end
            if n > 2
                d3Wdxi3 = dW{i}(:,:,3);
                dR{1}(:,counter,3) = (B{1}(:,k1,4).*W.^2 - 3*B{1}(:,k1,3).*W.^2.*dWdxi + B{1}(:,k1,2).*(-3*d2Wdxi2.*W.^2+6*W.*dWdxi.^2) ...
                                        + B{1}(:,k1,1).*(-d3Wdxi3.*W.^2 + 6*d2Wdxi2.*dWdxi.*W - 6*dWdxi.^3))./W.^4;
            end
            counter = counter + 1;
        end
    case 2
        for k2 = 1:degree(2)+1
            for k1 = 1:degree(1)+1    
                weight = weights(counter);

                W = W + B{1}(:,k1,1).*B{2}(:,k2,1)*weight;
                if n > 0
                    for i = 1:d_p
                        dW{i}(:,:,1) = dW{i}(:,:,1) + B{1}(:,k1,dd(i,1)).*B{2}(:,k2,dd(i,2))*weight;
                    end
                end
                counter = counter + 1;
            end
        end

        counter = 1;
        for k2 = 1:degree(2)+1
            for k1 = 1:degree(1)+1     
                fact = weights(counter)./(W.*W);

                BB = B{1}(:,k1,1).*B{2}(:,k2,1);
                R(:,counter) = BB.*fact.*W;

                if n > 0
                    for i = 1:d_p
                        dR{i}(:,counter,1) = (B{1}(:,k1,dd(i,1)).*B{2}(:,k2,dd(i,2)).*W - BB.*dW{i}(:,:,1)).*fact;
                    end
                end
                counter = counter + 1;
            end
        end
    case 3
        for k3 = 1:degree(3)+1
            for k2 = 1:degree(2)+1
                for k1 = 1:degree(1)+1   
                    weight = weights(counter);

                    W = W + B{1}(:,k1,1).*B{2}(:,k2,1).*B{3}(:,k3,1)*weight;
                    if n > 0
                        for i = 1:d_p
                            dW{i}(:,:,1) = dW{i}(:,:,1) + B{1}(:,k1,dd(i,1)).*B{2}(:,k2,dd(i,2)).*B{3}(:,k3,dd(i,3))*weight;
                        end
                    end
                    counter = counter + 1;
                end
            end
        end

        counter = 1;
        for k3 = 1:degree(3)+1
            for k2 = 1:degree(2)+1
                for k1 = 1:degree(1)+1       
                    fact = weights(counter)./(W.*W);

                    BB = B{1}(:,k1,1).*B{2}(:,k2,1).*B{3}(:,k3,1);
                    R(:,counter) = BB.*fact.*W;
                    if n > 0
                        for i = 1:d_p
                            dR{i}(:,counter,1) = (B{1}(:,k1,dd(i,1)).*B{2}(:,k2,dd(i,2)).*B{3}(:,k3,dd(i,3)).*W - BB.*dW{i}(:,:,1)).*fact;
                        end
                    end
                    counter = counter + 1;
                end
            end
        end
end
vargout = [{R},dR];