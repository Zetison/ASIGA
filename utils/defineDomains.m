function task = defineDomains(task)
task.varCol = copySet(task.varCol,1, 'inner', 'Gamma');
task.varCol = copySet(task.varCol,1, 'outer', 'Gamma_a');
task.varCol = copySet(task.varCol,1, 'inner', 'Neumann');
if numel(task.varCol) > 1
    task.varCol = copySet(task.varCol,2, 'outer', 'Gamma');
    task.varCol = copySet(task.varCol,2, 'outer', 'Neumann');
end
domain.Attributes.name = 'Omega_a';
domain.item{1}.Text = num2str(1:numel(task.varCol{1}.nurbs));
task.varCol{1}.geometry.domains.domain{1} = domain;