function [P_new, Xi_new] = insertKnotsInBsplines(n, p, Xi, newKnots, P)
% Implementation of Böhms method

Xi_new = Xi;
if isempty(newKnots)
    P_new = P;
    return;
else   
    P_new = zeros(size(P,1), n+length(newKnots));
end
for k = 1:length(newKnots)
    xi = newKnots(k);
	mu = findKnotSpan(n,p,xi,Xi_new);
    for i = 1:n+1
        if 1 <= i && i <= mu-p
            P_new(:,i) = P(:,i);
        elseif mu-p+1 <= i && i <= mu
            P_new(:,i) = (xi-Xi_new(i))/(Xi_new(i+p)-Xi_new(i))*P(:,i) ...
                 + (Xi_new(i+p)-xi)/(Xi_new(i+p)-Xi_new(i))*P(:,i-1);
        else
            P_new(:,i) = P(:,i-1);
        end
    end
    P = P_new(:,1:n+1);
    n = n+1;
    Xi_new = [Xi_new(1:mu) xi Xi_new(mu+1:end)];
end