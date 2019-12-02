function [Y, X] = difftable(A)
% 
%   
%    This is a solution of Hermite interpolation problem.
%   
%   example:
% 
%   A=[-1 2 -1 inf; 0 1 inf inf;1 -1 -1 8]
% 
%        
%     x    f(x)  f'(x)  f''(x)   . . .
%  
%  A =
%     -1     2    -1     Inf
%      0     1    Inf    Inf
%      1    -1    -1     8
% 
%   If you don't know the derive values, just write Inf.  
% 
%   Use this command: difftable(A)
%   And you can see the divided difference table, and the 
%   symbolic form of approximation polinomial.
% 
%   The function returns the coefficient vector of polinomial.
% 
%   You can use this function for calculate Newton form of 
%   interpolation. 
%    
%   example:
% 
%   A =
%     -1     2  
%      0     1  
%      1    -1  
% 
% 
%   Author:     Árpád Tóth
%               Eötvös University, Budapest
%               arpi@elte.hu

alappontok = A(:,1);
fvertekek = A(:,2);
[n,m] = size(A);
M = zeros(1,m*n,class(A));
M(1:n) = m-1;

k = sum(M);
B = zeros(k,k+1,class(A));           
C = zeros(k,size(A,2),class(A));           

l = 1;
aszam = 1;

for i = 1:k   
    if l > M(aszam)         
        l = 1;
        aszam = aszam + 1;       
    end
    B(i,1) = alappontok(aszam);
    B(i,2) = fvertekek(aszam); 
    C(i,:) = A(aszam,:);
    l = l + 1;
end

lvl = 2;
for l = 3:k+1
    I = 1:k+2-l; 
    indices = B(I,1) == B(I+l-2,1);
    II1 = find(indices);
    if ~isempty(II1)
        B(II1,l) = C(II1,l)/factorial(l-2);
    end
    II2 = find(~indices);
    if ~isempty(II2)
        B(II2,l) = (B(II2,lvl)-B(II2+1,lvl))./(B(II2,1)-B(II2+l-2,1));
    end
	lvl = lvl + 1;
end

Y = B(1,2:end);
X = B(:,1);


















