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
if strcmp(task.sol.solver,'LU')
    dA = decomposition(task.Pinv*task.A*task.Pinv,'lu');
end
for i = 1:noRHSs
    j = shiftROM + i-1; % Derivative order
    b = task.FF(:,i);   % Derivative of right hand side
    for k = 1:min(j,numel(dAdomega))
        b = b - nchoosek(j,k)*dAdomega{k}*task.U_p(:,shiftROM+i-k);
    end
    switch task.sol.solver
        case 'LU'
            task.U_p(:,shiftROM+i) = task.Pinv*(dA\(task.Pinv*b));
        case 'gmres'
            task.UU = zeros(size(task.FF));
            if task.sol.restart > size(task.A,1) % RESTART should be bounded by SIZE(task.A,1).
                task.sol.restart = size(task.A,1);
            end

            [task.U_p(:,shiftROM+i),task.sol.flag,task.sol.relres,task.sol.iter] = gmres(task.A,b,task.sol.restart,task.sol.tol,task.sol.maxit,task.L_A,task.U_A);
%                     gmres_parf_obj = parfeval(backgroundPool,@gmres,1,task.A,task.FF(:,i),task.sol.restart,task.sol.tol,task.sol.maxit,task.L_A,task.U_A);
%                     task.UU(:,i) = fetchOutputs(gmres_parf_obj);
        otherwise
            eval(['task.U_p(:,shiftROM+i) = ' task.sol.solver '(task.A,b,task.sol.tol,task.sol.maxit,task.L_A,task.U_A);'])
    end
end

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
