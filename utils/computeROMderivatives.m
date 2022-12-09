function task = computeROMderivatives(task,shiftROM)
% Construct the subspace V (scaled version of U = V/factorial(j))
if nargin < 2
    shiftROM = 0;
end
if task.misc.printLog
    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Computing derivatives for ROM ... ')
end

omega = task.misc.omega;
noRHSs = task.noRHSs;
dAdomega = cell(1,4);
dAdomega{1} = task.A1 + 2*omega*task.A2 + 4*omega^3*task.A4;
dAdomega{2} = 2*task.A2 + 12*omega^2*task.A4;
dAdomega{3} = 24*omega*task.A4;
dAdomega{4} = 24*task.A4;
if shiftROM > 0
    task.U = [task.U, zeros(size(task.FF))];
else
    task.U = zeros(size(task.FF));
end
dA = decomposition(task.Pinv*task.A*task.Pinv,'lu');
for i = 1:noRHSs
    j = shiftROM + i-1; % Derivative order
    b = task.FF(:,i);   % Derivative of right hand side
    for k = 1:min(j,numel(dAdomega))
        b = b - nchoosek(j,k)*dAdomega{k}*task.U(:,shiftROM+i-k);
    end
    task.U(:,shiftROM+i) = task.Pinv*(dA\(task.Pinv*b));
end

i = 1:noRHSs;
V = task.P_right*task.U(:,shiftROM+i);
task.V = [task.V, V./max(abs(V),[],1)];

if task.misc.printLog
    fprintf('using %12f seconds.', toc)
end