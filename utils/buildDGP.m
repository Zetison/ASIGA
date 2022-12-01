function task = buildDGP(task,startCol)
if nargin < 2
    startCol = 1;
end
if false
    % Perform modified Gram-Schmidt orthonormalization
    k = size(task.V,2);
    for i = startCol:k
        task.V(:,i) = task.V(:,i) / norm(task.V(:,i));
        for j = i+1:k
            task.V(:,j) = task.V(:,j) - proj(task.V(:,i),task.V(:,j));
        end
    end
else
    for i = startCol:size(task.V,2)
        V_i = task.V(:,i);
        for j = 1:i-1
            task.V(:,i) = task.V(:,i) - dot(task.V(:,j),V_i)/dot(task.V(:,j),task.V(:,j))*task.V(:,j);
        end
        task.V(:,i) = task.V(:,i)/norm(task.V(:,i));
    end
end

% Build reduced order matrices
task.A0_am = task.V'*task.A0*task.V;
task.A1_am = task.V'*task.A1*task.V;
task.A2_am = task.V'*task.A2*task.V;
task.A4_am = task.V'*task.A4*task.V;

function w = proj(u,v)
w = (dot(v,u) / dot(u,u)) * u;