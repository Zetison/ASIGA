function A = invVandermonde(x)
% Formula at https://proofwiki.org/wiki/Inverse_of_Vandermonde%27s_Matrix

n = length(x);
% A = zeros(n);
% 
% for i = 1:n
%     for j = 1:n
%         denominator = 1;
%         for m = setdiff(1:n, j)
%             denominator = denominator*(x(m) - x(j));
%         end
%         m_values = nchoosek(setdiff(1:n, j),n-i);
%         numerator = 0;
%         for k = size(m_values,1)
%             numerator = numerator + (-1)^i*prod(x(m_values(k,:)));
%         end
%         A(i,j) = numerator/denominator;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LU decomposition given in http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19660023042.pdf
L_inv = zeros(n);
U_inv = zeros(n);

L_inv(1,1) = 1;
for i = 2:n
    for j = 1:i
        L_inv(i,j) = 1/prod(x(j) - x(setdiff(1:i, j)));
    end
end

for i = 1:n
    U_inv(i,i) = 1;
end
for j = 2:n
    U_inv(1,j) = -U_inv(1,j-1)*x(j-1);
end
for i = 2:n
    for j = i+1:n
        U_inv(i,j) = U_inv(i-1,j-1) - U_inv(i,j-1)*x(j-1);
    end
end
            
A = U_inv*L_inv;