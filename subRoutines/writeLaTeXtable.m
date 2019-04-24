coreMethods = {'IGA','IGA','hp_FEM','linear_FEM'};
BC = {'SHBC','SSBC','NNBC'};

for i_coreM = 1:3
    fid = fopen(['../articles/articleIGA/contents/IGAFEM_' BC{i_coreM} '.tex'], 'w+','b');
    meshStrings = {'${\\cal M}_{6,1,\\mathrm{i}}^{\\textsc{fem}}$', ...
                   '${\\cal M}_{5,2,\\mathrm{i}}^{\\textsc{fem}}$', ...
                   '${\\cal M}_{5,2,1}^{\\textsc{iga}}$\t\t', ...
                   '${\\cal M}_{4,3,2}^{\\textsc{iga}}$\t\t', ...
                   '${\\cal M}_{5,3,2}^{\\textsc{iga}}$\t\t'};
    fprintf(fid, ['Mesh ' meshStrings{1} '\t\t&']);
    j = 9+i_coreM;
    jj = 2;
    printRow(fid,studies,j,jj)
    
    fprintf(fid, ['Mesh ' meshStrings{2} '\t\t&']);
    j = 6+i_coreM;
    jj = 2;
    printRow(fid,studies,j,jj)
    
    fprintf(fid, ['Mesh ' meshStrings{3} '\t\t&']);
    j = 3+i_coreM;
    jj = 3;
    printRow(fid,studies,j,jj)
    
    fprintf(fid, ['Mesh ' meshStrings{4} '\t\t&']);
    j = i_coreM;
    jj = 2;
    printRow(fid,studies,j,jj)
    
    fprintf(fid, ['Mesh ' meshStrings{5} '\t\t&']);
    j = 3+i_coreM;
    jj = 4;
    printRow(fid,studies,j,jj)
    fclose(fid); 
end

function printRow(fid,studies,j,jj)
fprintf(fid, ' %d\t& %d\t& %3.2f\t& %3.2f\t& %3.2f\t\\cr\n', studies(j).results(jj).varCol.noElems, ...
                                                             studies(j).results(jj).varCol.dofs, ...
                                                             studies(j).results(jj).varCol.timeBuildSystem, ...
                                                             studies(j).results(jj).varCol.tot_time-studies(j).results(jj).varCol.timeBuildSystem, ...
                                                             studies(j).results(jj).energyError);
end