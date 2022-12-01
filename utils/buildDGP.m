function task = buildDGP(task,startCol)
if nargin < 2
    startCol = 1;
end
for i = startCol:size(task.V,2)
    for j = 1:i-1
        task.V(:,i) = task.V(:,i) - dot(task.V(:,j),task.V(:,i))*task.V(:,j);
    end
    task.V(:,i) = task.V(:,i)/norm(task.V(:,i));
end

% Build reduced order matrices
task.A0_am = task.V'*task.A0*task.V;
task.A1_am = task.V'*task.A1*task.V;
task.A2_am = task.V'*task.A2*task.V;
task.A4_am = task.V'*task.A4*task.V;

function w = proj(u,v)
w = (dot(v,u) / dot(u,u)) * u;