function printResultsToTextFiles(study,newOptions)

options = struct('xname',           'alpha',  ...
                 'yname',           'TS', ...
                 'plotResults', 	1, ... 
                 'printResults',	0, ... 
                 'axisType',        'plot', ... 
                 'lineStyle',        '*-', ... 
                 'xScale',          1, ...
                 'yScale',          1, ... 
                 'xlim',            NaN, ...
                 'ylim',            NaN, ...
                 'xLoopName',       NaN, ...
                 'legendEntries',   {{}}, ...
                 'subFolderName',   '', ...
                 'fileDataHeaderX', [], ...
                 'noXLoopPrms',     0);
if nargin > 1
    options = updateOptions(options,newOptions);
end
xname = options.xname;
yname = options.yname;
plotResults = options.plotResults;
printResults = options.printResults;
legendEntries = options.legendEntries;
if ~plotResults && ~printResults
    return
end
loopParametersArr = study.loopParametersArr;
loopParameters = study.loopParameters;
noParms = numel(loopParameters);
model = study.tasks(1).task.model;
sizes = zeros(1,length(loopParametersArr));
for i = 1:length(loopParametersArr)
    sizes(i) = numel(loopParametersArr{i});
end
% sizes = fliplr(sizes);
noTasks = prod(sizes);
if isempty(options.subFolderName)
    subFolderName = study.resultsFolder;
else
    subFolderName = options.subFolderName;
end
if ~exist(subFolderName, 'dir')
    mkdir(subFolderName);
end
analyticSolutionExist = study.tasks(1).task.varCol{1}.analyticSolutionExist;
% analyticSolutionExist = true;

if isempty(options.fileDataHeaderX)
    fileDataHeaderX = xname;
else
    fileDataHeaderX = options.fileDataHeaderX;
end
if ~isempty(get(gca, 'Children'))
    legendPrev = legend;
    legendArr = legendPrev.String;
    analyticPlotted = false;
    for i = 1:numel(legendArr)
        if strcmp(legendArr{i},'Analytic solution')
            analyticPlotted = true;
        end
    end
else
    analyticPlotted = false;
end

plotAnalyticSolution = analyticSolutionExist && (strcmp(yname, 'TS') || strcmp(yname, 'p') || strcmp(yname, 'abs_p')) && ~analyticPlotted;
switch options.noXLoopPrms
    case 0
        if plotAnalyticSolution
            y_ref = study.tasks(1).task.results.([yname '_ref']);
            x = study.tasks(1).task.varCol{1}.(xname);
            x = x*options.xScale;
            y_ref = y_ref*options.yScale;
            plotXY(x,y_ref,options.axisType,options.lineStyle,[0,0,0],'Analytic solution');
            hold on
        end
        col = LaTeXcolorMap(noTasks);
        for i = 1:noTasks
            if isfield(study.tasks(i).task.varCol{1},xname)
                x = study.tasks(i).task.varCol{1}.(xname);
            else
                x = study.tasks(i).task.(xname);
            end
            x = x*options.xScale;
            if ~isfield(study.tasks(i).task.results, yname)
                y = NaN*x';
            else
                y = study.tasks(i).task.results.(yname);
            end
            y = y*options.yScale;
            
            [legendName, saveName] = constructStrings(legendEntries,i,i,study,model,noParms,loopParameters,yname,fileDataHeaderX);
            if printResults
                if isrow(y)
                    y = y.';
                end
                saveName(saveName == '.') = [];
                printResultsToFile2([subFolderName '/' saveName], {'x', x.', 'y', y, 'xlabel',{fileDataHeaderX}, 'ylabel',{yname}, ...
                                                                                                            'task', study.tasks(i).task});
            end
            if plotResults
%                 legendName(legendName == '_') = [];
                plotXY(x,y,options.axisType,options.lineStyle,col(i,:),legendName);
                hold on
            end
        end
        if plotAnalyticSolution && printResults
            if isrow(y_ref)
                y_ref = y_ref.';
            end
            printResultsToFile2([subFolderName '/' saveName], {'x', x.', 'y', y_ref, 'xlabel',{fileDataHeaderX}, 'ylabel',{yname}, ...
                                                                                                            'task', study.tasks(i).task});
        end
    case 1
        idx = 1;
        while ~strcmp(loopParameters{idx},options.xLoopName)
            idx = idx + 1;
            if idx > noParms
                error('xname is not one of the loopParameters')
            end
        end
        x = zeros(sizes); % x axis data
        y = zeros(sizes); % y axis data 
        if plotAnalyticSolution
            y_ref = zeros(sizes); % y axis data
        end
        idxMap = zeros(sizes); % y axis data 

        for i = 1:noTasks
            if isfield(study.tasks(i).task.varCol{1},xname)
                x(i) = study.tasks(i).task.varCol{1}.(xname);  
            else
                x(i) = study.tasks(i).task.(xname);  
            end  
            if ~isfield(study.tasks(i).task.results, yname)
                y(i) = NaN*x(i);
            else
                y(i) = study.tasks(i).task.results.(yname);
            end
            if plotAnalyticSolution
                y_ref(i) = study.tasks(i).task.results.([yname '_ref']);
            end
            idxMap(i) = i;
        end
        % Collect error data
        otherInd = [1:idx-1,idx+1:noParms];
        otherParms = sizes(otherInd);
        
        x = permute(x, [setdiff(1:noParms,idx), idx])*options.xScale;
        y = permute(y, [setdiff(1:noParms,idx), idx])*options.yScale;
        if plotAnalyticSolution
            y_ref = permute(y_ref, [setdiff(1:noParms,idx), idx])*options.yScale;
        end

        idxMap = permute(idxMap, [setdiff(1:noParms,idx), idx]);

        size1 = sizes(idx);
        size2 = prod(otherParms);
        col = LaTeXcolorMap(size2);
        for i = 1:size2
            indices = i:size2:(size1-1)*size2+i;
            x_temp = x(indices);
            y_temp = y(indices);
            if plotAnalyticSolution
                y_ref_temp = y_ref(indices);
            end

            [legendName, saveName] = constructStrings(legendEntries,i,idxMap,study,model,noParms,loopParameters,yname,fileDataHeaderX,otherInd);
            if printResults
                saveName(saveName == '.') = [];
                printResultsToFile2([subFolderName '/' saveName], {'x', x_temp(:), 'y', y_temp(:), 'xlabel',{fileDataHeaderX}, 'ylabel',{yname}, ...
                                                                                'task', study.tasks(idxMap(i)).task, 'xLoopName',options.xLoopName});
            end
            if plotResults
                x_temp = x_temp(:);
                y_temp = y_temp(:);
                if plotAnalyticSolution
                    y_ref_temp = y_ref_temp(:);
                    plotXY(x_temp,y_ref_temp,options.axisType,options.lineStyle,[0,0,0],'Analytic solution');
                    hold on
                end
%                 legendName(legendName == '_') = [];
                plotXY(x_temp,y_temp,options.axisType,options.lineStyle,col(i,:),legendName);
                hold on
            end
        end
        if plotAnalyticSolution && printResults
            y_ref = permute(y_ref,[noParms,1:noParms-1]);
            y_ref = y_ref(1:sizes(idx));
            x = permute(x,[noParms,1:noParms-1]);
            x = x(1:sizes(idx));
            if isrow(y_ref)
                y_ref = y_ref.';
            end
            printResultsToFile2([subFolderName '/' saveName], {'x', x.', 'y', y_ref, 'xlabel',{fileDataHeaderX}, 'ylabel',{yname}, ...
                                                                'task', study.tasks(i).task});
        end
end

if plotResults
    title(['Results for model ' study.tasks(1).task.model],'interpreter','none')
    intrprtrX = 'latex';
    intrprtrY = 'latex';
    switch xname
        case 'alpha'
            xLabel = '$$\alpha$$';
        case {'omega','omega_ROM'}
            xLabel = '$$\omega$$';
        case {'f','f_ROM'}
            xLabel = '$$f$$';
        case {'k','k_ROM'}
            xLabel = '$$k$$';
        case 'beta'
            xLabel = '$$\beta$$';
        case 's_ie'
            xLabel = '$$s$$';
        case 'dofs'
            xLabel = 'Degrees of freedom';
        case {'surfDofs'}
            xLabel = 'Dofs at $$\Gamma$$';
        case 'dofsAlg'
            xLabel = 'Degrees of freedom in algebraic scale $$N^{1/3}$$';
            xt = get(gca, 'XTick');
            set(gca, 'XTick', xt, 'XTickLabel', xt.^3)
        case 'nepw'
            xLabel = 'Number of elements per wave length';
        case 'tot_time'
            xLabel = 'Total time';
        case 'timeSolveSystem'
            xLabel = 'Time solving the system';
        case 'timeBuildSystem'
            xLabel = 'Time building the system';
        case 'N'
            xLabel = '$$N$$';
        case 'agpBEM'
            xLabel = '$$s$$';
        otherwise
            xLabel = xname;
            intrprtrX = 'none';
    end
    if (strcmp(xname,'alpha') || strcmp(xname,'beta')) && (options.xScale == 180/pi)
        if strcmp(options.axisType,'polar')
            thetatickformat('degrees')
        else
            xtickformat('degrees')
        end
    end
    switch yname
        case 'p'
            yLabel = '$$p$$';
        case 'TS'
            yLabel = '$$\mathrm{TS}$$';
        case 'abs_p'
            yLabel = '$$|p|$$';
        case 'error_p'
            yLabel = '$$\frac{|p-p_h|}{|p|}$$';
        case 'error_pAbs'
            yLabel = '$$\frac{||p|-|p_h||}{|p|}$$';
        case 'surfaceError'
            yLabel = '$$\frac{\|p-p_h\|_{L_2(\Gamma)}}{\|p\|_{L_2(\Gamma)}}$$';
        case 'energyError'
            yLabel = '$$\frac{\|p-p_h\|_E}{\|p\|_E}$$';
        case 'H1sError'
            yLabel = '$$\frac{|p-p_h|_{H^1}}{|p|_{H^1}}$$';
        case 'L2Error'
            yLabel = '$$\frac{\|p-p_h\|_{L^2}}{\|p\|_{L^2}}$$';
        case 'cond_number'
            yLabel = 'Condition number';
        otherwise
            yLabel = yname;
            intrprtrY = 'none';
    end
    if (strcmp(yname,'error_pAbs') || strcmp(yname,'error_p') || strcmp(yname,'surfaceError') || strcmp(yname,'energyError')) && (options.yScale == 1)
        yLabel = [yLabel, ' [\%]'];
    end
            
          
    if ~strcmp(options.axisType,'polar')  
        ylabel(yLabel,'interpreter',intrprtrY)
        xlabel(xLabel,'interpreter',intrprtrX)
    end
    h = gca;
    noGraphs = numel(h.Children);
    col = LaTeXcolorMap(noGraphs);
    for i = 1:noGraphs
        h.Children(i).Color = col(i,:);
    end
    if ~isnan(options.xlim)
        if strcmp(options.axisType,'polar')  
            thetalim(options.xlim)
        else
            xlim(options.xlim)
        end
    end
    if ~isnan(options.ylim)
        if strcmp(options.axisType,'polar')  
            rlim(options.ylim)
        else
            ylim(options.ylim)
        end
    end
    savefig([subFolderName '/plot_' model '_' yname 'VS' xname])
end

function [legendName, saveName] = constructStrings(legendEntries,i,idxMap,study,model,noParms,loopParameters,yname,xname,otherInd)
saveName = model;           
legendName = [];          
            
if isempty(legendEntries)
    if nargin > 9
        for j = otherInd
            temp2 = loopParameters{j};
            temp = study.tasks(idxMap(i)).task.(loopParameters{j});
            if isnumeric(temp) || islogical(temp)
                temp = num2str(temp);
            end
            if isstruct(temp)
                fieldNames = fieldnames(temp);
                temp2 = fieldNames{1};
                temp = temp.(temp2);
            end
            [mathematicalLegend, legendEntriesMath, scale, postfix] = fixLegendEntry(loopParameters{j});
            [saveName, legendName] = updateStrings(saveName, legendName, j, loopParameters, mathematicalLegend, legendEntriesMath, scale, postfix,temp,temp2);
        end
    else
        for j = 1:noParms
            temp2 = loopParameters{j};

            temp = study.tasks(i).task.(loopParameters{j});
            if ~ischar(temp)
                temp = num2str(temp);
            end
            if isstruct(temp)
                fieldNames = fieldnames(temp);
                temp2 = fieldNames{1};
                temp = temp.(temp2);
            end
            [mathematicalLegend, legendEntriesMath, scale, postfix] = fixLegendEntry(loopParameters{j});
            [saveName, legendName] = updateStrings(saveName, legendName, j, loopParameters, mathematicalLegend, legendEntriesMath, scale, postfix,temp,temp2);
        end 
    end  
else
    for j = 1:length(legendEntries)
        temp2 = legendEntries{j};
        if ~isfield(study.tasks(i).task,legendEntries{j})
            continue
        end
        temp = study.tasks(i).task.(legendEntries{j});
        if strcmp(temp,'NaN')
            continue
        end
        [mathematicalLegend, legendEntriesMath, scale, postfix] = fixLegendEntry(legendEntries{j});
        [saveName, legendName] = updateStrings(saveName, legendName, j, legendEntries, mathematicalLegend, legendEntriesMath, scale, postfix,temp,temp2);
    end
end
saveName = [saveName '_' yname 'VS' xname];

function [saveName, legendName] = updateStrings(saveName, legendName, j, legendEntries, mathematicalLegend, legendEntriesMath, scale, postfix,temp,temp2)
if isnumeric(temp) || islogical(temp)
    temp = num2str(temp*scale);
end
if isstruct(temp)
    fieldNames = fieldnames(temp);
    temp2 = fieldNames{1};
    temp = temp.(temp2);
end
if ~isempty(legendName)
    legendName = [legendName ', ' ];
end
switch legendEntries{j}
    case {'formulation','IEbasis','method','coreMethod','BC'}
        legendName = [legendName temp];
        saveName = [saveName '_' temp];
    otherwise
        saveName = [saveName '_' temp2 temp];
        if mathematicalLegend
            legendName = [legendName '$$' legendEntriesMath '$$=' temp postfix];
        else
            legendEntry = legendEntries{j};
            legendEntry = insertBefore(legendEntry,'_','\');
            legendName = [legendName legendEntry '=' temp];
        end
end

function [mathematicalLegend, legendEntriesMath, scale, postfix] = fixLegendEntry(legendEntry)
mathematicalLegend = false;
scale = 1;
postfix = '';
legendEntriesMath = NaN;
switch legendEntry
    case 'extraGP'
        mathematicalLegend = true;
        legendEntriesMath = 'n_{\mathrm{eq}}^{(1)}';
    case 'extraGPBEM'
        mathematicalLegend = true;
        legendEntriesMath = 'n_{\mathrm{eq}}^{(2)}';
    case 'colBEM_C0'
        mathematicalLegend = true;
        legendEntriesMath = 'C_{\mathrm{col}}';
    case 'c_x'
        mathematicalLegend = true;
        legendEntriesMath = 'c_x';
    case 'agpBEM'
        mathematicalLegend = true;
        legendEntriesMath = 's';
    case 's_ie'
        mathematicalLegend = true;
        legendEntriesMath = 's';
    case 'alpha_s'
        mathematicalLegend = true;
        legendEntriesMath = '{\alpha}_s';
        postfix = '$$^\circ$$';
        scale = 180/pi;
    case 'beta_s'
        mathematicalLegend = true;
        legendEntriesMath = '\beta';
        postfix = '$$^\circ$$';
        scale = 180/pi;
end