function task = buildDGP(task,startCol)
if nargin < 2
    startCol = 1;
end
% if task.rom.adaptiveROM
%     if 1
%         if false
%             [~, indices] = sort(task.rom.J_U);
%             [~, indices] = sort(task.rom.omega_U(indices));
%         else
%             [~, indices] = sort(task.rom.omega_U);
%         end
%         task.V = task.U(:,indices);
%     else
%         task.V = task.U;
%     end
% else
% %     [~, indices] = sort(task.rom.omega_U);
% %     task.V = task.U(:,indices);
%     task.V = task.U;
% end
algorithm = 1;
switch algorithm
    case 1
        % Perform Gram-Schmidt orthogonalization
%         for i = 1:size(task.V,2)
        for i = startCol:size(task.V,2)
            for j = 1:i-1
                task.V(:,i) = task.V(:,i) - (task.V(:,j)'*task.V(:,i))*task.V(:,j);
            end
            task.V(:,i) = task.V(:,i)/norm(task.V(:,i));
        end
    case 2
        % https://blogs.mathworks.com/cleve/2016/07/25/compare-gram-schmidt-and-householder-orthogonalization-algorithms/
        U = task.V;
        [n,p] = size(U);
        task.V = zeros(n,p);
        R = zeros(p,p);
        for k = 1:p
            task.V(:,k) = U(:,k);
            for i = 1:k-1
                R(i,k) = task.V(:,i)'*task.V(:,k);
                task.V(:,k) = task.V(:,k) - R(i,k)*task.V(:,i);
            end
            R(k,k) = norm(task.V(:,k))';
            task.V(:,k) = task.V(:,k)/R(k,k);
        end
    case 3
        % https://github.com/areslp/matlab/blob/master/drtoolbox/techniques/mgs.m
        U = task.V;
        [~, n] = size(U);
        V = U;
        R = zeros(n, n);
        for i=1:n
            R(i,i) = norm(V(:,i));
            V(:,i) = V(:,i) / R(i, i);
            if (i < n)
                for j=i+1:n
                    R(i,j) = V(:,i)' * V(:,j);
                    V(:,j) = V(:,j) - R(i, j) * V(:,i);
                end
            end
        end
        task.V = V;
    case 4
        % https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
        U = task.V;
        [n,k] = size(U);
        task.V = zeros(n,k);
        task.V(:,1) = U(:,1)/norm(U(:,1));
        for i = 2:k
            task.V(:,i) = U(:,i);
            for j=1:i-1
                task.V(:,i)=task.V(:,i)-(task.V(:,j)'*task.V(:,i)) * task.V(:,j);
            end
            task.V(:,i) = task.V(:,i)/norm(task.V(:,i));
        end
end

% Build reduced order matrices
task.A0_am = task.V'*(task.P_rightinv'*task.A0*task.P_rightinv)*task.V;
task.A1_am = task.V'*(task.P_rightinv'*task.A1*task.P_rightinv)*task.V;
task.A2_am = task.V'*(task.P_rightinv'*task.A2*task.P_rightinv)*task.V;
if task.misc.symmetric
    task.A4_am = task.V'*(task.P_rightinv'*task.A4*task.P_rightinv)*task.V;
else
    task.A4_am = zeros(size(task.A0_am));
end

