function task = ASIGAparaview(task)
resultsFolder = task.resultsFolder;
printLog = task.misc.printLog;
if ~task.rom.useROM && ~strcmp(task.misc.method,'RT')
    tic
    if task.misc.plotResidualError && task.analyticSolutionExist && strcmp(task.misc.formulation,'GCBIE')
%                 close all
        fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Plotting zeros of residual ... ')
        plotGalerkinResidual(task.varCol{1});
        fprintf('using %12f seconds.', toc)
        axis equal
        axis off
%             set(gca, 'Color', 'none');
%                 title('Error in Galerkin solution')
%                 colorbar
%                 caxis([-7,0])
%                 camlight
        colorbar off
        savefig([resultsFolder '/' task.saveName])
%             n_xi = fluid.number(1);
%             n_eta = fluid.number(2);
%             [cg_xi, grev_xi] = CauchyGalerkin(fluid.degree(1), n_xi, fluid.knots{1});
%             [cg_eta, grev_eta] = CauchyGalerkin(fluid.degree(2), n_eta, fluid.knots{2});
%             % keyboard
%             n_cp = task.varCol{1}.noDofs - length(dofsToRemove);
%             
%             cp_cg = zeros(n_cp,3);
%             cp_grev = zeros(n_cp,3);
%             counter = 1;
%             counter2 = 1;
%             for j = 1:n_eta
%                 eta_cg = cg_eta(j);
%                 eta_grev = grev_eta(j);
%                 for i = 1:n_xi
%                     if ~any(dofsToRemove == counter)
%                         xi_cg = cg_xi(i);        
%                         xi_grev = grev_xi(i);            
%                         cp_cg(counter2,:) = evaluateNURBS(fluid, [xi_cg, eta_cg]); 
%                         cp_grev(counter2,:) = evaluateNURBS(fluid, [xi_grev, eta_grev]); 
%                         counter2 = counter2 + 1;
%                     end
%                     counter = counter + 1;
%                 end
%             end
%             hold on
%             h1 = plot3(cp_cg(:,1),cp_cg(:,2),cp_cg(:,3), '*','color','red');
%             h2 = plot3(cp_grev(:,1),cp_grev(:,2),cp_grev(:,3), '*','color','blue');
%             legend([h1, h2], {'CG','Grev'})
%             hold off
%             savefig([resultsFolderName '/' task.saveName '_surfPlot_mesh' num2str(M) '_formulation_' formulation '_degree' num2str(max(fluid.degree)) '.fig'])
    end
    para = task.para;
    if para.plotResultsInParaview
        task.ffp.alpha_s = task.ffp.alpha_s(task.para.i_MS);
        task = getAnalyticSolutions(task);
        if isempty(para.name)
            resultsFolderNameParaview = [resultsFolder '/paraviewResults'];
            if ~exist(resultsFolderNameParaview, 'dir')
                mkdir(resultsFolderNameParaview);
            end          
            vtfFileName = [resultsFolderNameParaview '/' task.saveName];
            para.name = vtfFileName;
        end
        M = task.msh.M;

        para.extraXiPts = eval(para.extraXiPts);
        para.extraEtaPts = eval(para.extraEtaPts);
        para.extraZetaPts = eval(para.extraZetaPts);
        if strcmp(task.misc.method, 'KDT')
            task.varCol{1}.U = zeros(task.varCol{1}.noDofs,1);
        end
        if para.plotFullDomain
            tic
            if printLog
                fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Creating volumetric paraview files for all domains ... ')
            end
            createParaviewFiles(task, 'para_options', para);
        end

        for i_v = 1:numel(task.varCol)
            if ~isempty(para.plotSubsets) && ~task.varCol{i_v}.boundaryMethod
                name = para.name;
                for i = 1:numel(para.plotSubsets)
                    bdryName = para.plotSubsets{i};
                    if printLog
                        fprintf(['\n%-' num2str(task.misc.stringShift) 's'], ['Creating paraview files for ' bdryName ' in domain ' num2str(i_v) ' ... '])
                    end
                    [~,setFound] = findSet(task.varCol{i_v}.geometry.topologysets.set,bdryName);
                    if setFound
                        taskBdry = task;
                        taskBdry.varCol{i_v} = rmfield(taskBdry.varCol{i_v},'geometry');
                        varColBdry = meshBoundary(task.varCol{i_v},bdryName);
                        taskBdry.varCol{i_v}.nurbs = varColBdry.nurbs;
                        taskBdry = collectVariables(taskBdry);
                        taskBdry.varCol{i_v} = findDofsToRemove(generateIGAmesh(taskBdry.varCol{i_v}));
    %                     fieldNames = fields(varColBdry);
    %                     for j = 1:numel(fieldNames)
    %                         taskBdry.varCol{i_v}.(fieldNames{j}) = varColBdry.(fieldNames{j});
    %                     end
                        nodes = varColBdry.nodes;

                        d_f = task.varCol{i_v}.dimension;
                        nodesField = zeros(1,d_f*numel(nodes));
                        for j = 1:d_f
                            nodesField(j:d_f:end) = d_f*(nodes-1)+j;
                        end
                        taskBdry.varCol{i_v}.U = task.varCol{i_v}.U(nodesField,:);

                        para.name = [name, '_', bdryName];

                        createParaviewFiles(taskBdry, 'para_options', para);
                    else
                        switch bdryName
                            case {'xy','yz','xz'}
                                warning(['The set ' bdryName ' was not found for domain ' num2str(i_v) ', and is thus not plotted in paraview. For xy, yz and xz; setting msh.explodeNURBS = true might help.'])
                            otherwise
                                if ~(strcmp(bdryName,'Gamma_a') && i_v > 1)
                                    warning(['The set ' bdryName ' was not found for domain ' num2str(i_v) ', and is thus not plotted in paraview.'])
                                end
                        end
                    end
                end
            end
        end
    end
end