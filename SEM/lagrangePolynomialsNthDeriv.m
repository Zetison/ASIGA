function dp = lagrangePolynomialsNthDeriv(x,indices,N,X,D,d)
if nargin < 6
    if isrow(x)
        x = x';
    end
    d = 1;
end

dp = zeros(numel(x),D-d+1);
for l = 1:N
    if ~ismember(l,indices)
        j = setdiff(1:N,[indices,l]);        
        dp(:,1) = dp(:,1) + prod((x-X(j))./repmat(X(indices(1))-X(j),numel(x),1),2)/(X(indices(1))-X(l));
        if d < D
            dp(:,2:end) = dp(:,2:end) + lagrangePolynomialsNthDeriv(x,[indices,l],N,X,D,d+1)/(X(indices(1))-X(l));
        end
    end
end