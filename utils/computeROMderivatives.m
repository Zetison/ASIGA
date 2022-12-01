function task = computeROMderivatives(task,shiftROM)
if nargin < 2
    shiftROM = 0;
end
if task.misc.printLog
    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Computing derivatives for ROM ... ')
end

omega = task.misc.omega;
dAdomega = cell(1,4);
dAdomega{1} = task.A1 + 2*omega*task.A2 + 4*omega^3*task.A4;
dAdomega{2} = 2*task.A2 + 12*omega^2*task.A4;
dAdomega{3} = 24*omega*task.A4;
dAdomega{4} = 24*task.A4;
noColsInV = size(task.V,2);
task.V = [task.V, zeros(size(task.FF))];
dA = decomposition(task.Pinv*task.A*task.Pinv,'lu');
for i = 1:task.noRHSs
    j = shiftROM + i-1;
    b = task.FF(:,i);
    for k = 1:min(j,numel(dAdomega))
        b = b - nchoosek(j,k)*dAdomega{k}*task.V(:,noColsInV+i-k);
    end
    task.V(:,noColsInV+i) = task.Pinv*(dA\(task.Pinv*b));
end

if task.misc.printLog
    fprintf('using %12f seconds.', toc)
end