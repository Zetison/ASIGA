function plotROMresiduals(task)
c_f = task.varCol{1}.c_f;
if isfield(task.rom,'history')
    history = task.rom.history;
    for i = 1:numel(history)
        figure(20+i)
        omega = task.rom.history(i).omega;
        residual = task.rom.history(i).residual;
        hold on
    
        omega_P = history(i).omega_P;
        J_P = history(i).J_P;
        semilogy(omega_P/c_f,1e-15*ones(size(omega_P)),'o','color','magenta','DisplayName','Interpolation points')
        for j = 1:numel(omega_P)
            text((omega_P(j)+omega_P(end)/200)/c_f,1e-15,num2str(J_P(j)),'color','magenta')
        end

        semilogy([omega_P(1),omega_P(end)]/c_f,task.rom.tolerance*ones(1,2),'red','DisplayName','Tolerance')

        semilogy(omega/c_f,residual,'*','color','black','DisplayName','Residual check points')

        omega_T_new = history(i).omega_T_new;
        [~,idx] = ismember(omega,omega_T_new);
        semilogy(omega_T_new/c_f,residual(logical(idx)),'o','color','cyan','DisplayName','New interpolation point')

        ylim([5e-16,200])
        set(gca,'yscale','log')
        legend show
        savefig([task.resultsFolder, '/', task.saveName '_adaptiveROM_iter' num2str(i)])
    end
    legend off
    if isfield(task.rom.history(end),'residualFine')
        residual = task.rom.history(end).residualFine;
        omegaFine = task.rom.history(end).omegaFine;
        semilogy(omegaFine/c_f,residual,'blue','DisplayName','Final relative residual')
    end
    if isfield(task.rom.history(end),'relError')
        relError = task.rom.history(end).relError;
        semilogy(omegaFine/c_f,relError,'green','DisplayName','Final relative error')
    end
    set(gca,'yscale','log')
    ylim([5e-16,200])
    ylabel('Relative error/residual')
    xlabel('Wavenumber')
    legend show
    savefig([task.resultsFolder, '/', task.saveName '_adaptiveROM_iter' num2str(i)])
end