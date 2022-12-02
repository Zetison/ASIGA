function task = Hetmaniuk2012raa_P1(task)


if task.misc.printLog
    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Computing ROM basis adaptively ... ')
end
omega = task.misc.omega;
n_c = task.rom.n_c;
oldprintLog = task.misc.printLog;
task.misc.printLog = false;

% Algorithm P1 in Hetmaniuk2013aas available at https://www.doi.org/10.1002/nme.4436
task.V = [];
omega_T = [];
omega_T_new = [task.rom.omega(1),task.rom.omega(end)];
task.rom.history = struct();
J_P = [];
counter = 1;
while ~isempty(omega_T_new)
    J_new = zeros(1,numel(omega_T_new));
    for p = 1:numel(omega_T_new)
        % Algorithm P2 in Hetmaniuk2013aas
        [task, J_new(p)] = Hetmaniuk2012raa_P2(task, union(omega_T_new,omega_T), omega_T_new(p));
    end
    [omega_T,sortIdx] = sort([omega_T,omega_T_new]);
    omega_T_new = [];
    J_P = [J_P,J_new];
    J_P = J_P(sortIdx);
    task.rom.history(counter).residual = [];
    task.rom.history(counter).omega = [];
    for p = 1:numel(omega_T)-1
        omega_p = omega_T(p);
        omega_pp1 = omega_T(p+1);
        j = 1:n_c;
        omega_cj = omega_p + j/(n_c+1)*(omega_pp1 - omega_p);
        residual = computeROMresidual(task, omega_cj);
        task.rom.history(counter).residual = [task.rom.history(counter).residual, residual];
        task.rom.history(counter).omega = [task.rom.history(counter).omega, omega_cj];
        [~,j_max] = max(residual);
        if residual(j_max) > task.rom.tolerance
            omega_T_new = union(omega_T_new,omega_cj(j_max));
        end
    end
    task.rom.history(counter).omega_P = omega_T;
    task.rom.history(counter).J_P = J_P;
    counter = counter + 1;
end
if task.rom.computeROMresidualFine
    residual = computeROMresidual(task, omega);
    task.rom.history(end).residualFine = residual;
%     semilogy(omega, residual,'DisplayName',num2str(counter))
%     hold on
%     set(gca,'yscale','log')
end
% close all
% history = task.rom.history;
% for i = 1:numel(history)
%     semilogy(history(i).omega, history(i).residual,'DisplayName',num2str(i))
%     hold on
% end
% legend show
% set(gca,'yscale','log')
% keyboard
task.misc.printLog = oldprintLog;
if task.misc.printLog
    fprintf('using %12f seconds.', toc)
end