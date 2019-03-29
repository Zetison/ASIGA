function dp = lagrangePolynomialsNthDerivOld(x,indices,N,Xi,D,d)
if isrow(x)
    x = x';
end
if nargin < 6
    d = 1;
end
% if nargin < 7
%     dp = zeros(numel(x),D);
% end

dp = zeros(numel(x),1);
for l = 1:N
    if ~ismember(l,indices)
        j = setdiff(1:N,[indices,l]);
%         dpdxi = dpdxi + prod((xi-Xi(l))./repmat(Xi(i)-Xi(l),numel(xi),1),2)/(Xi(i)-Xi(j));
        
%         dp(:,d:end) = dp(:,d:end) + repmat(prod((x-Xi(j))./repmat(Xi(indices(1))-Xi(j),numel(x),1),2)/(Xi(indices(1))-Xi(l)),1,D-d+1);
%         if d < D
%             dp = dp + lagrangePolynomialsNthDeriv(x,[indices,l],N,Xi,D,d+1,dp)/(Xi(indices(1))-Xi(l));
%         end
        if d == D
            dp = dp + prod((x-Xi(j))./repmat(Xi(indices(1))-Xi(j),numel(x),1),2)/(Xi(indices(1))-Xi(l));
        else
            dp = dp + lagrangePolynomialsNthDerivOld(x,[indices,l],N,Xi,D,d+1)/(Xi(indices(1))-Xi(l));
        end
    end
end