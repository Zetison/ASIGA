function task = buildDGP(task,startCol)
if nargin < 2
    startCol = 1;
end

% Perform Gram-Schmidt orthogonalization
% for i = 1:size(task.V,2)
for i = startCol:size(task.V,2)
    for j = 1:i-1
        task.V(:,i) = task.V(:,i) - (task.V(:,j)'*task.V(:,i))*task.V(:,j);
    end
    task.V(:,i) = task.V(:,i)/norm(task.V(:,i));
end

% Build reduced order matrices
task.A0_am = task.V'*task.A0*task.V;
task.A1_am = task.V'*task.A1*task.V;
task.A2_am = task.V'*task.A2*task.V;
task.A4_am = task.V'*task.A4*task.V;


return


% https://github.com/areslp/matlab/blob/master/drtoolbox/techniques/mgs.m
% [m, n] = size(task.U);
% V = task.U;
% R = zeros(n, n);
% for i=1:n
%     R(i,i) = norm(V(:,i));
%     V(:,i) = V(:,i) / R(i, i);
%     if (i < n)
%         for j=i+1:n
%             R(i,j) = V(:,i)' * V(:,j);
%             V(:,j) = V(:,j) - R(i, j) * V(:,i);
%         end
%     end
% end
% task.V = V;

% % https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
% [n,k] = size(task.U);
% task.V = zeros(n,k);
% task.V(:,1) = task.U(:,1)/norm(task.U(:,1));
% for i = 2:k
%     task.V(:,i) = task.U(:,i);
%     for j=1:i-1
%         task.V(:,i)=task.V(:,i)-(task.V(:,j)'*task.V(:,i)) * task.V(:,j);
%     end
%     task.V(:,i) = task.V(:,i)/norm(task.V(:,i));
% end
% 
% https://blogs.mathworks.com/cleve/2016/07/25/compare-gram-schmidt-and-householder-orthogonalization-algorithms/
% [n,p] = size(task.U);
% task.V = zeros(n,p);
% R = zeros(p,p);
% for k = 1:p
%     task.V(:,k) = task.U(:,k);
%     for i = 1:k-1
%         R(i,k) = task.V(:,i)'*task.V(:,k);
%         task.V(:,k) = task.V(:,k) - R(i,k)*task.V(:,i);
%     end
%     R(k,k) = norm(task.V(:,k))';
%     task.V(:,k) = task.V(:,k)/R(k,k);
% end