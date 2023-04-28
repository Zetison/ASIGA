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
    task.U_p = [task.U_p, zeros(size(task.FF))];
else
    task.U_p = zeros(size(task.FF));
end
tic
if true
    dA = decomposition(task.Pinv*task.A*task.Pinv,'lu');
    for i = 1:noRHSs
        j = shiftROM + i-1; % Derivative order
        b = task.FF(:,i);   % Derivative of right hand side
        for k = 1:min(j,numel(dAdomega))
            b = b - nchoosek(j,k)*dAdomega{k}*task.U_p(:,shiftROM+i-k);
        end
        task.U_p(:,shiftROM+i) = task.Pinv*(dA\(task.Pinv*b));
    end
else
    dA = decomposition(task.Pinv*task.A*task.Pinv,'lu');
    for i = 1:noRHSs
        j = shiftROM + i-1; % Derivative order
        b = task.FF(:,i);   % Derivative of right hand side
        for k = 1:min(j,numel(dAdomega))
            b = b - nchoosek(j,k)*dAdomega{k}*task.U_p(:,shiftROM+i-k);
        end
        task.U_p(:,shiftROM+i) = task.Pinv*(dA\(task.Pinv*b));
    end
end
toc

i = 1:noRHSs;
% if ~isfield(task,'U')
%     task.U = [];
% end
% if 1
%     if 0
%         task.U = [task.U, task.U_p(:,shiftROM+i)];
%     else
%         V = task.P_right*task.U_p(:,shiftROM+i);
%         task.U = [task.U, V./max(abs(V),[],1)];
%     end
%     V = task.P_right*task.U_p(:,shiftROM+i);
%     task.V = [task.V, V./max(abs(V),[],1)];
% else
%     task.V = [task.V, task.U_p(:,shiftROM+i)];
%     task.U = [task.U, task.P_right*task.U_p(:,shiftROM+i)];
% end
V = task.P_right*task.U_p(:,shiftROM+i);
task.V = [task.V, V./max(abs(V),[],1)];

if task.misc.printLog
    fprintf('using %12f seconds.', toc)
end