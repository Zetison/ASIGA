function boundaryMethod = isBoundaryMethod(task)

switch task.misc.method
    case {'IE','ABC','PML'}
        boundaryMethod = false;
    case {'IENSG'}
        boundaryMethod = task.iem.boundaryMethod;
    case {'BEM','KDT','MFS','RT'}
        boundaryMethod = true;
    case 'BA'
        switch task.misc.formulation
            case {'SL2E','SL2Etot'}
                boundaryMethod = true;
            case {'VL2E','VL2Etot'}
                boundaryMethod = false;
            otherwise
                error('Not implemented')
        end
end