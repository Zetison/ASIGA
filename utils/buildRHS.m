function task = buildRHS(task,overrideUseROM)

if nargn < 2
    overrideUseROM = false;
end
if overrideUseROM
    oldUseROM = task.rom.useROM;
    oldNoRHSs = task.noRHSs;
    task.rom.useROM = true;
    task.noRHSs = 1;
end

noDomains = numel(task.varCol);
for i_domain = 1:min(noDomains,2)
    task = applyNeumannCondition(task,i_domain);
    if task.misc.symmetric && strcmp(task.varCol{i_domain}.media,'solid')
        task.varCol{i_domain}.FF = omega_i^2*task.varCol{i_domain}.FF;
    end
end
if overrideUseROM
    task.rom.useROM = oldUseROM;
    task.noRHSs = oldNoRHSs;
end
