function [P_new, Xi_new] = elevateBsplinesDegree(n, p, Xi, P_old, m)
% Implementation of method in "Efficient Degree Elevation and Knot Insertion
% for B-spline Curves using Derivatives" by Qi-Xing Huang a Shi-Min Hu, 
% Ralph R Martin.


P = zeros(size(P_old,1),n,p+1);
S = numel(unique(Xi)) - 1;
n_new = n + S*m;

% Find multiplicity vector Xi
z = ones(1,S+1);
l = 1;
for i = 1:n+p
    if Xi(i+1) == Xi(i)
        z(l) = z(l)+1;
    else
        l = l+1;
    end
end

    
% Find P_i^j
P(:,:,1) = P_old;
for j = 1:p
    for i = 1:n
        if Xi(i+p+1) > Xi(i+j)
            P(:,i,j+1) = (p+1-j)/(Xi(i+p+1)-Xi(i+j))*(P(:,i+1,j) - P(:,i,j));
        end
    end
end

% Find new knot vector
Xi_new = zeros(1,n_new+p+m+1);
temp = unique(Xi);
j = 1;
for i = 1:length(temp)
    Xi_new(j:j+z(i)-1+m) = temp(i)*ones(1,z(i)+m);
    j = j+z(i)+m;
end
Q = zeros(size(P,1), n_new, p+1);
for j = 0:p
    Q(:,1,j+1) = P(:,1,j+1);
end
for i = 1:S-1
    beta = 0;
    for l = 1:i
        beta = beta + z(l+1);
    end
    for j = p+1-z(i+1):p
        Q(:,beta+1+i*m,j+1) = P(:,beta+1,j+1);
    end
end
for i = 0:S-1
    beta = 0;
    for l = 1:i
        beta = beta + z(l+1);
    end
    for k = 1:m
        Q(:,beta+1+i*m+k,p+1) = Q(:,beta+1+i*m,j+1);
    end
%     for k = m+1:m-1+z(i+1+1)
%         Q(:,beta+1+i*m+k,p+1) = zeros(size(Q,1),1);
%     end
end
for j = p:-1:1
    for i = 1:n_new-1
        if Xi_new(i+p+m+1) > Xi_new(i+j)
            Q(:,i+1,j-1+1) = Q(:,i,j-1+1) + (Xi_new(i+1+p+m)-Xi_new(i+j))/(p+m+1-j)*Q(:,i,j+1);
        end
    end
end

P_new = Q(:,:,1);
        