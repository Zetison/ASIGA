function task = buildDGP(task)

% Perform Gram-Schmidt orthonormalization
task.V(:,1) = task.V(:,1)/sqrt(task.V(:,1)'*task.V(:,1));
for i = 2:size(task.V,2)
    task.V(:,i) = task.V(:,i);
    for j = 1:i-1
        task.V(:,i) = task.V(:,i) - ( task.V(:,j)'*task.V(:,i) )/( task.V(:,j)'*task.V(:,j) )*task.V(:,j);
    end
    task.V(:,i) = task.V(:,i)/sqrt(task.V(:,i)'*task.V(:,i));
end

% Build reduced order matrices
task.A0_am = task.V'*task.A0*task.V;
task.A1_am = task.V'*task.A1*task.V;
task.A2_am = task.V'*task.A2*task.V;
task.A4_am = task.V'*task.A4*task.V;