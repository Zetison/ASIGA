function [U,dU] = projectBC(varCol,SHBC,useCBIE,useHBIE)

dofsToRemove = varCol.dofsToRemove;
if SHBC
    U = NaN;
    dU = NaN;
    if useCBIE
        varCol2 = varCol;
        varCol2.formulation = 'SL2E';
        varCol2.analytic = varCol.p_inc;
        [A, FF] = bestApproximationVec(varCol2);
        noDofs_tot = varCol.noDofs;

        A(dofsToRemove,:) = [];  
        A(:,dofsToRemove) = [];

        FF(dofsToRemove,:) = [];
        UU = A\FF;
        U = zeros(noDofs_tot,size(UU,2),size(UU,3));
        U(setdiff(1:noDofs_tot, dofsToRemove'),:,:) = UU; 
    end
    if useHBIE
        varCol2 = varCol;
        varCol2.formulation = 'SL2E';
        varCol2.analytic = varCol.dp_inc;
        [A, FF] = bestApproximationVec(varCol2);
        noDofs_tot = varCol.noDofs;

        A(dofsToRemove,:) = [];  
        A(:,dofsToRemove) = [];

        FF(dofsToRemove,:) = [];
        dUU = A\FF;
        dU = zeros(noDofs_tot,size(dUU,2),size(dUU,3));
        dU(setdiff(1:noDofs_tot, dofsToRemove'),:,:) = dUU; 
    end
else
    varCol2 = varCol;
    varCol2.formulation = 'SL2E';
    varCol2.analytic = varCol.dpdn;
    [A, FF] = bestApproximationVec(varCol2);
    noDofs_tot = varCol.noDofs;

    A(dofsToRemove,:) = [];  
    A(:,dofsToRemove) = [];

    FF(dofsToRemove,:) = [];
    UU = A\FF;
    U = zeros(noDofs_tot,size(UU,2),size(UU,3));
    U(setdiff(1:noDofs_tot, dofsToRemove'),:,:) = UU;   
    dU = NaN;
end