function printResultsToTextFiles(study,newOptions)

options = struct('xname',           'alpha',  ...
                 'yname',           'TS', ...
                 'plotResults', 	1, ... 
                 'printResults',	0, ... 
                 'addSlopes',	    0, ... 
                 'axisType',        'plot', ... 
                 'lineStyle',        '*-', ... 
                 'xScale',          1, ...
                 'yScale',          1, ... 
                 'xLabel',          [], ...
                 'yLabel',          [], ...
                 'format',          'LaTeX', ...
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
model = study.tasks(1).task.misc.model;
sizes = zeros(1,length(loopParametersArr));
for i = 1:length(loopParametersArr)
    sizes(i) = numel(loopParametersArr{i});
end
if length(loopParametersArr) == 1
    sizes = [sizes,1];
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

if isempty(options.fileDataHeaderX)
    fileDataHeaderX = xname;
else
    fileDataHeaderX = options.fileDataHeaderX;
end
ioOptions.format = options.format;

switch options.noXLoopPrms
    case 0
        col = LaTeXcolorMap(noTasks);
        for i = 1:noTasks
            try 
                x = eval(['study.tasks(i).task.' xname]);
            catch
                x = NaN;
                warning('ASIGA:xDataNotFound', [xname ' was not found in the structure.'])
            end
            x = x*options.xScale;
            try 
                y = eval(['study.tasks(i).task.results.' yname]);
            catch
                y = NaN*x';
                warning('ASIGA:yDataNotFound', [yname ' was not found in the structure.'])
            end
            y = y*options.yScale;
            
            [legendName, saveName] = constructStrings(legendEntries,i,i,study,model,noParms,loopParameters,yname,fileDataHeaderX);
            if printResults
                if isrow(y)
                    y = y.';
                end
                if isrow(x)
                    x = x.';
                end
                saveName(saveName == '.') = [];
                saveName(saveName == '\') = [];
                saveName(saveName == '/') = [];
                
                ioOptions.filename = [subFolderName '/' saveName];
                ioOptions.x = x;
                ioOptions.y = y;
                ioOptions.xlabel = {fileDataHeaderX};
                ioOptions.ylabel = {yname};
                ioOptions.task = study.tasks(i).task;
                ioOptions.xLoopName = NaN;
                printResultsToFile2(ioOptions);
            end
            analyticSolutionExist = study.tasks(i).task.analyticSolutionExist;
            plotAnalyticSolution = analyticSolutionExist && (strcmp(yname, 'TS') || strcmp(yname, 'p') || strcmp(yname, 'p_Re') || ...
                                                                                strcmp(yname, 'p_Im') || strcmp(yname, 'abs_p'));
            if plotAnalyticSolution
                y_ref = study.tasks(i).task.results.([yname '_ref']);
                y_ref = y_ref*options.yScale;
            end
            if plotResults
                if plotAnalyticSolution
                    prevGrafs = get(gca, 'Children');
                    analyticPlotted = false;
                    if ~isempty(prevGrafs)
                        for j = 1:numel(prevGrafs)
                            if numel(y_ref) == numel(prevGrafs(j).YData) && all(y_ref(:) == prevGrafs(j).YData(:))
                                analyticPlotted = true;
                            end
                        end
                    end
                    if ~analyticPlotted
                        y_ref = study.tasks(i).task.results.([yname '_ref']);
                        y_ref = y_ref*options.yScale;
                        plotXY(x,y_ref,options.axisType,options.lineStyle,[0,0,0],'Analytic solution',options.addSlopes);
                        hold on
                    end
                end
                
                plotXY(x,y,options.axisType,options.lineStyle,col(i,:),legendName,options.addSlopes);
                hold on
            end
            if plotAnalyticSolution && printResults
                if isrow(y_ref)
                    y_ref = y_ref.';
                end
                if isrow(x)
                    x = x.';
                end
                ioOptions.filename = [subFolderName '/' saveName '_analytic'];
                ioOptions.x = x;
                ioOptions.y = y_ref;
                ioOptions.xlabel = {fileDataHeaderX};
                ioOptions.ylabel = {yname};
                ioOptions.task = study.tasks(i).task;
                ioOptions.xLoopName = NaN;
                printResultsToFile2(ioOptions);
            end
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
        analyticSolutionExist = study.tasks(1).task.analyticSolutionExist;
        plotAnalyticSolution = analyticSolutionExist && (strcmp(yname, 'TS') || strcmp(yname, 'p') || strcmp(yname, 'p_Re') || strcmp(yname, 'p_Im') || strcmp(yname, 'abs_p'));
        if plotAnalyticSolution
            y_ref = zeros(sizes); % y axis data
        end
        idxMap = zeros(sizes); % y axis data 

        for i = 1:noTasks
            try 
                x(i) = eval(['study.tasks(i).task.' xname]);
            catch
                warning('ASIGA:xDataNotFound', [xname ' was not found in the structure.'])
            end
            try 
                y(i) = eval(['study.tasks(i).task.results.' yname]);
            catch
                y(i) = NaN;
                warning('ASIGA:yDataNotFound', [yname ' was not found in the structure.'])
            end
            if plotAnalyticSolution
                y_ref(i) = study.tasks(i).task.results.([yname '_ref']);
            end
            idxMap(i) = i;
        end
        % Collect error data
        otherInd = [1:idx-1,idx+1:noParms];
        otherParms = sizes(otherInd);
        if noParms > 1
            x = permute(x, [setdiff(1:noParms,idx), idx])*options.xScale;
            y = permute(y, [setdiff(1:noParms,idx), idx])*options.yScale;
            if plotAnalyticSolution
                y_ref = permute(y_ref, [setdiff(1:noParms,idx), idx])*options.yScale;
            end
            idxMap = permute(idxMap, [setdiff(1:noParms,idx), idx]);
        end
        
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
                saveName(saveName == '\') = [];
                saveName(saveName == '/') = [];
                
                ioOptions.filename = [subFolderName '/' saveName];
                ioOptions.x = x_temp(:);
                ioOptions.y = y_temp(:);
                ioOptions.xlabel = {fileDataHeaderX};
                ioOptions.ylabel = {yname};
                ioOptions.task = study.tasks(idxMap(i)).task;
                ioOptions.xLoopName = options.xLoopName;
                printResultsToFile2(ioOptions);
            end
            if plotResults
                x_temp = x_temp(:);
                y_temp = y_temp(:);
                if plotAnalyticSolution
                    prevGrafs = get(gca, 'Children');
                    analyticPlotted = false;
                    if ~isempty(prevGrafs)
                        for j = 1:numel(prevGrafs)
                            if numel(y_ref_temp) == numel(prevGrafs(j).YData) && all(y_ref_temp(:) == prevGrafs(j).YData(:))
                                analyticPlotted = true;
                            end
                        end
                    end
                    if ~analyticPlotted
                        y_ref_temp = y_ref_temp(:);
                        plotXY(x_temp,y_ref_temp,options.axisType,options.lineStyle,[0,0,0],'Analytic solution',options.addSlopes);
                        hold on
                    end
                end
                plotXY(x_temp,y_temp,options.axisType,options.lineStyle,col(i,:),legendName,options.addSlopes);
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
            if isrow(x)
                x = x.';
            end
                
            ioOptions.filename = [subFolderName '/' saveName '_analytic'];
            ioOptions.x = x;
            ioOptions.y = y_ref;
            ioOptions.xlabel = {fileDataHeaderX};
            ioOptions.ylabel = {yname};
            ioOptions.task = study.tasks(i).task;
            printResultsToFile2(ioOptions);
        end
end

if plotResults
    title(['Results for model ' study.tasks(1).task.misc.model],'interpreter','none')
    intrprtrX = 'latex';
    intrprtrY = 'latex';
    if isempty(options.xLabel)
        switch xname
            case 'ffp.alpha'
                xLabel = '$$\alpha$$';
            case {'misc.omega'}
                xLabel = '$$\omega$$';
            case {'misc.f'}
                xLabel = '$$f$$';
            case {'varCol{1}.k'}
                if options.xScale == 1
                    xLabel = '$$k$$';
                else
                    xLabel = '$$ka$$';
                end
            case 'ffp.beta'
                xLabel = '$$\beta$$';
            case 'iem.s_ie'
                xLabel = '$$s$$';
            case 'dofs'
                xLabel = 'Degrees of freedom';
            case 'surfDofs'
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
            case 'iem.N'
                xLabel = '$$N$$';
            case 'pml.gamma'
                xLabel = '$$\gamma$$';
            case 'bem.agpBEM'
                xLabel = '$$s$$';
            case 'varCol{1}.kL'
                xLabel = '$$kL$$';
            case 'varCol{1}.ka'
                xLabel = '$$ka$$';
            case 'varCol{1}.kR'
                xLabel = '$$kR$$';
            otherwise
                xLabel = xname;
                intrprtrX = 'none';
        end
    else
        xLabel = options.xLabel;
    end
    if (strcmp(xname,'ffp.alpha') || strcmp(xname,'ffp.beta')) && (options.xScale == 180/pi)
        if strcmp(options.axisType,'polar')
            thetatickformat('degrees')
        else
            xtickformat('degrees')
        end
    end
    if isempty(options.yLabel)
        switch yname
            case 'p'
                yLabel = '$$p$$';
            case 'p_Re'
                yLabel = 'Real part of pressure';
            case 'p_Im'
                yLabel = 'Imaginary part of pressure';
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
    else
        yLabel = options.yLabel;
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
    idx = find(xname == '.',1,'last');
    if ~isempty(idx)
        xname = xname(idx+1:end);
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
            if numel(idxMap) > 1
                temp = eval(['study.tasks(idxMap(i)).task.' loopParameters{j}]);
            else
                temp = eval(['study.tasks(idxMap).task.' loopParameters{j}]);
            end
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

            temp = eval(['study.tasks(i).task.' loopParameters{j}]);
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
        try
            if numel(idxMap) > 1
                temp = eval(['study.tasks(idxMap(i)).task.' legendEntries{j}]);
            else
                temp = eval(['study.tasks(idxMap).task.' legendEntries{j}]);
            end
            if strcmp(temp,'NaN')
                continue
            end
            [mathematicalLegend, legendEntriesMath, scale, postfix] = fixLegendEntry(legendEntries{j});
            [saveName, legendName] = updateStrings(saveName, legendName, j, legendEntries, mathematicalLegend, legendEntriesMath, scale, postfix,temp,temp2);
        catch ME
            error([legendEntries{j} ' is not a field of task'])
        end
    end
end
idx = find(xname == '.',1,'last');
if isempty(idx)
    saveName = [saveName '_' yname 'VS' xname];  
else  
    saveName = [saveName '_' yname 'VS' xname(idx+1:end)];
end

function [saveName, legendName] = updateStrings(saveName, legendName, j, legendEntries, mathematicalLegend, legendEntriesMath, scale, postfix,temp,temp2)
if isnumeric(temp) || islogical(temp)
    temp = num2str(temp*scale);
end
if isstruct(temp)
    fieldNames = fieldnames(temp);
    temp2 = fieldNames{1};
    temp = temp.(temp2);
end
idx = find(temp2 == '.',1,'last');
if ~isempty(idx)
    temp2 = temp2(idx+1:end);
end
if ~isempty(legendName)
    legendName = [legendName ', ' ];
end
switch legendEntries{j}
    case {'misc.scatteringCase','misc.formulation','ie.IEbasis','misc.method','misc.coreMethod','misc.BC'}
        temp = insertBefore(temp,'_','\');
        legendName = [legendName temp];
        saveName = [saveName '_' temp];
    otherwise
        saveName = [saveName '_' temp2 temp];
        if mathematicalLegend
            legendName = [legendName '$$' legendEntriesMath '$$=' temp postfix];
        else
            legendEntry = legendEntries{j};
            legendEntry = insertBefore(legendEntry,'_','\');
            idx = find(legendEntry == '.',1,'last');
            if ~isempty(idx)
                legendEntry = legendEntry(idx+1:end);
            end
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
    case 'gamma'
        mathematicalLegend = true;
        legendEntriesMath = '\gamma';
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