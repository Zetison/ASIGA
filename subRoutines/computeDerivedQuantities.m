function task = computeDerivedQuantities(task)


totNoElems = 0;
FEdofs = 0;
for i = 1:task.noDomains
    totNoElems = totNoElems + task.varCol{i}.noElems;
    FEdofs = FEdofs + task.varCol{i}.noDofs - numel(task.varCol{i}.dofsToRemove);
end
task.totNoElems = totNoElems;
task.FEdofs = FEdofs;
if task.misc.compute_h_max
    h_max = -Inf;
    for i = 1:task.noDomains
        if task.msh.refineThetaOnly
            h_max_i = findMaxElementDiameter(subNURBS(task.varCol{i}.nurbs,'at',[1,0;0,0;0,0]));
        else
            h_max_i = findMaxElementDiameter(task.varCol{i}.nurbs);
        end
        if h_max < h_max_i
            h_max = h_max_i;
        end
    end
    task.h_max = h_max;
    if isfield(task.varCol{1},'c_f')
        task.varCol{1}.nepw = task.varCol{1}.lambda./task.h_max;
    end
end

task.surfDofs = getNoSurfDofs(task);
for field = {'a','R','L'}
    if isfield(task.varCol{1},field{1})
        task.varCol{1}.(['k' field{1}]) = task.varCol{1}.(field{1})*task.varCol{1}.k;
    end
end
if isfield(task.misc,'omega')
    task.f = task.misc.omega/(2*pi);
end
if isfield(task,'dofs')
    task.dofsAlg = (task.dofs)^(1/3);
end