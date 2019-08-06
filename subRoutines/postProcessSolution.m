
%% Sort data from U into U_fluid_o, U_solid, and U_fluid_i
if strcmp(method,'MFS')
    U_fluid_o = UU;
else
    U = zeros(noDofs_tot,size(UU,2),size(UU,3));
    U(setdiff(1:noDofs_tot, dofsToRemove'),:,:) = UU;   

    if ~useSolidDomain
    %         U_fluid_o = U(1:varCol_fluid_o.noDofs,:,:); 
        U_fluid_o = U; 
        U_fluid_o = addSolutionToRemovedNodes_new(U_fluid_o, varCol);
    elseif useSolidDomain && ~useInnerFluidDomain
        U_fluid_o = U((varCol_solid.noDofs+1):end,:,:);
        U_solid = U(1:varCol_solid.noDofs,:,:); 

        U_fluid_o = addSolutionToRemovedNodes_new(U_fluid_o, varCol);
        U_solid = addSolutionToRemovedNodes_new(U_solid, varCol_solid);
    else
        U_fluid_o = U((varCol_fluid_i.noDofs+varCol_solid.noDofs+1):end,:,:); 
        U_solid = U((varCol_fluid_i.noDofs+1):(varCol_fluid_i.noDofs + varCol_solid.noDofs),:,:);
        U_fluid_i = U(1:varCol_fluid_i.noDofs,:,:);

        U_fluid_o = addSolutionToRemovedNodes_new(U_fluid_o, varCol);
        U_solid = addSolutionToRemovedNodes_new(U_solid, varCol_solid);
        U_fluid_i = addSolutionToRemovedNodes_new(U_fluid_i, varCol_fluid_i);
    end
    % save([resultsFolderName '/' saveName  '_' coreMethod method '_mesh' num2str(M) '_degree' num2str(max(varCol.nurbs.degree)) ...
    %                         '_formulation' varCol.formulation '.mat'], 'varCol', 'U_fluid_o')
end

%% Find h_max and store results
if i_k == 1
    h_max_fluid = findMaxElementDiameter(varCol.patches);
    h_max = h_max_fluid;
    if useSolidDomain
        h_max_solid = findMaxElementDiameter(varCol_solid.patches);
        h_max = max([h_max, h_max_solid]);
    end
    if useInnerFluidDomain
        h_max_fluid_i = findMaxElementDiameter(varCol_fluid_i.patches);
        h_max = max([h_max, h_max_fluid_i]);
    end
    varCol.h_max = h_max;
    varCol.nepw = lambda(1)./h_max;
    varCol.dofs = actualNoDofs;
    varCol.surfDofs = getNoSurfDofs(varCol);
    if storeSolution
        varCol.U = U_fluid_o;
    end
end