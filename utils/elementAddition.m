function P = elementAddition(v, e)
% v is a 1xM matrix and e is a NxM matrx. Computes the sums 
% P(i,:) = v + e(i,:).
% N = size(e,1);
% P = zeros(size(e));
% for i = 1:N
%     P(i,:) = v + e(i,:);
% end
P = repmat(v,size(e,1),1)+e;