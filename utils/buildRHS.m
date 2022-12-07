function task = buildRHS(task,overrideUseROM,noRHSs)

if nargin < 2
    overrideUseROM = false;
end
if nargin < 3
    noRHSs = 1;
end
if overrideUseROM
    oldUseROM = task.rom.useROM;
    oldNoRHSs = task.noRHSs;
    task.rom.useROM = false;
    task.noRHSs = noRHSs;
end

noDomains = numel(task.varCol);
for i_domain = 1:min(noDomains,2)
    task = applyNeumannCondition(task,i_domain);
end
if overrideUseROM
    task.rom.useROM = oldUseROM;
    task.noRHSs = oldNoRHSs;
end
