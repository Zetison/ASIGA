function printResultsToTextFiles(study,newOptions)

options = struct('xname',           'alpha',  ...
                 'yname',           'TS', ...
                 'plotResults', 	1, ... 
                 'printResults',	0, ... 
                 'axisType',        'plot', ... 
                 'lineStyle',        '*-', ... 
                 'xScale',          1, ...
                 'yScale',          1, ... 
                 'xLoopName',       NaN, ...
                 'legendEntries',   {{}}, ...
                 'subFolderName',   NaN, ...
                 'fileDataHeaderX', [], ...
                 'noXLoopPrms',     0);
             
if nargin > 1
    newOptionFields = fieldnames(newOptions);
    for j = 1:numel(newOptionFields)
        options.(newOptionFields{j}) = newOptions.(newOptionFields{j});
    end
end
xname = options.xname;
yname = options.yname;
plotResults = options.plotResults;
printResults = options.printResults;
legendEntries = options.legendEntries;


if ~exist(study.subFolderName, 'dir')
    mkdir(study.subFolderName);
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
if isnan(options.subFolderName)
    subFolderName = study.subFolderName;
else
    subFolderName = options.subFolderName;
end
analyticSolutionExist = study.tasks(1).task.varCol.analyticSolutionExist;
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
            x = study.tasks(1).task.varCol.(xname);
            x = x*options.xScale;
            y_ref = y_ref*options.yScale;
            plotXY(x,y_ref,options.axisType,options.lineStyle,[0,0,0],'Analytic solution');
            hold on
        end
        col = LaTeXcolorMap(noTasks);
        for i = 1:noTasks
            if isfield(study.tasks(i).task.varCol,xname)
                x = study.tasks(i).task.varCol.(xname);
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
            
            saveName = model;           
            legendName = [];
            for j = 1:noParms
                temp2 = loopParameters{j};

                temp = study.tasks(i).task.(loopParameters{j});
                if ~ischar(temp)
                    temp = num2str(temp);
                end
                if j ~= 1
                    legendName = [legendName ', '];
                end
                if isstruct(temp)
                    fieldNames = fieldnames(temp);
                    temp2 = fieldNames{1};
                    temp = temp.(temp2);
                end
                switch loopParameters{j}
                    case {'formulation','IEbasis','method','coreMethod','BC'}
                        saveName = [saveName '_' temp];
                        legendName = [legendName temp];
                    otherwise
                        saveName = [saveName '_' temp2 temp];
                        legendName = [legendName loopParameters{j} '=' temp];
                end
            end 
            if ~isempty(legendEntries)
                legendName = [];
                saveName = model;  
                for j = 1:length(legendEntries)
                    temp2 = legendEntries{j};
                    if ~isfield(study.tasks(i).task,legendEntries{j})
                        continue
                    end
                    temp = study.tasks(i).task.(legendEntries{j});
                    if isnumeric(temp)
                        temp = num2str(temp);
                    end
                    if j ~= 1
                        legendName = [legendName ', ' ];
                    end
                    if isstruct(temp)
                        fieldNames = fieldnames(temp);
                        temp2 = fieldNames{1};
                        temp = temp.(temp2);
                    end
                    switch legendEntries{j}
                        case {'formulation','IEbasis','method','coreMethod','BC'}
                            legendName = [legendName temp];
                            saveName = [saveName '_' temp];
                        otherwise
                            saveName = [saveName '_' temp2 temp];
                            legendName = [legendName legendEntries{j} '=' temp];
                    end
                end
            end
            if printResults
                if isrow(y)
                    y = y.';
                end
                saveName(saveName == '.') = [];
                printResultsToFile2([subFolderName '/' saveName], x.', y, {fileDataHeaderX}, {yname}, study.tasks(i).task);
            end
            if plotResults
                legendName(legendName == '_') = [];
                plotXY(x,y,options.axisType,options.lineStyle,col(i,:),legendName);
                xlim([min(x),max(x)])
                hold on
            end
        end
        if plotAnalyticSolution && printResults
            if isrow(y_ref)
                y_ref = y_ref.';
            end
            printResultsToFile2([subFolderName '/analytic'], x.', y_ref, {fileDataHeaderX}, {yname}, study.tasks(i).task);
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
            if isfield(study.tasks(i).task.varCol,xname)
                x(i) = study.tasks(i).task.varCol.(xname);  
            else
                x(i) = study.tasks(i).task.(xname);  
            end  
            y(i) = study.tasks(i).task.results.(yname);
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

            saveName = model;
            legendName = [];

            for j = otherInd
                temp2 = loopParameters{j};
                temp = study.tasks(idxMap(i)).task.(loopParameters{j});
                if isnumeric(temp)
                    temp = num2str(temp);
                end
                if j ~= otherInd(1)
                    legendName = [legendName ', ' ];
                end
                if isstruct(temp)
                    fieldNames = fieldnames(temp);
                    temp2 = fieldNames{1};
                    temp = temp.(temp2);
                end
                switch loopParameters{j}
                    case {'formulation','IEbasis','method','coreMethod','BC'}
                        saveName = [saveName '_' temp];
                        legendName = [legendName temp];
                    otherwise
                        saveName = [saveName '_' temp2 temp];
                        legendName = [legendName loopParameters{j} '=' temp];
                end
            end
            if ~isempty(legendEntries)
                legendName = [];
                saveName = model;  
                for j = 1:length(legendEntries)
                    temp2 = legendEntries{j};
                    if ~isfield(study.tasks(idxMap(i)).task,legendEntries{j})
                        continue
                    end
                    temp = study.tasks(idxMap(i)).task.(legendEntries{j});
                    if isnumeric(temp)
                        temp = num2str(temp);
                    end
                    mathematicalLegend = false;
                    switch legendEntries{j}
                        case 'extraGP'
                            mathematicalLegend = true;
                            legendEntriesMath = 'n_{\mathrm{eq}}^{(1)}';
                        case 'extraGPBEM'
                            mathematicalLegend = true;
                            legendEntriesMath = 'n_{\mathrm{eq}}^{(2)}';
                        case 'colBEM_C0'
                            mathematicalLegend = true;
                            legendEntriesMath = 'C_{\mathrm{col}}';
                        case 'agpBEM'
                            mathematicalLegend = true;
                            legendEntriesMath = 's';
                    end
                    if isstruct(temp)
                        fieldNames = fieldnames(temp);
                        temp2 = fieldNames{1};
                        temp = temp.(temp2);
                    end
                    if strcmp(temp,'NaN')
                        continue
                    end
                    if j ~= 1
                        legendName = [legendName ', ' ];
                    end
                    switch legendEntries{j}
                        case {'formulation','IEbasis','method','coreMethod','BC'}
                            legendName = [legendName temp];
                            saveName = [saveName '_' temp];
                        otherwise
                            saveName = [saveName '_' temp2 temp];
                            if mathematicalLegend
                                legendName = [legendName '$' legendEntriesMath '$=' temp];
                            else
                                legendName = [legendName legendEntries{j} '=' temp];
                            end
                    end
                end
            end
            if printResults
                saveName(saveName == '.') = [];
                printResultsToFile2([subFolderName '/' saveName], x_temp(:), y_temp(:), {fileDataHeaderX}, {yname}, study.tasks(idxMap(i)).task, options.xLoopName);
            end
            if plotResults
                x_temp = x_temp(:);
                y_temp = y_temp(:);
                xlimOld = xlim;
                if xlimOld(1) == 0 && xlimOld(2) == 1
                    xlimOld = [inf,-inf];
                end
                if plotAnalyticSolution
                    y_ref_temp = y_ref_temp(:);
                    plotXY(x_temp,y_ref_temp,options.axisType,options.lineStyle,[0,0,0],'Analytic solution');
                    hold on
                end
                legendName(legendName == '_') = [];
                plotXY(x_temp,y_temp,options.axisType,options.lineStyle,col(i,:),legendName);
                xLim = [min(min(x_temp),xlimOld(1)),max(max(x_temp),xlimOld(2))];
                if xLim(1) ~= xLim(2)
                    xlim(xLim)
                end
                hold on
            end
        end
        if plotAnalyticSolution && printResults
            y_ref = permute(y_ref,[noParms,1:noParms-1]);
            y_ref = y_ref(1:sizes(idx));
            x = permute(x,[noParms,1:noParms-1]);
            x = x(1:sizes(idx));
            printResultsToFile2([subFolderName '/analytic'], x, y_ref, {fileDataHeaderX}, {yname}, study.tasks(i).task);
        end
end

if plotResults
    title(['Results for model ' study.tasks(1).task.model],'interpreter','none')
    intrprtrX = 'latex';
    intrprtrY = 'latex';
    switch xname
        case 'alpha'
            xLabel = '$$\alpha$$';
        case 'f'
            xLabel = '$$f$$';
        case {'k','k_ROM'}
            xLabel = '$$k$$';
        case 'beta'
            xLabel = '$$\beta$$';
        case 'dofs'
            xLabel = 'Degrees of freedom';
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
        xtickformat('degrees')
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
            yLabel = 'Conditioning number';
        otherwise
            yLabel = yname;
            intrprtrY = 'none';
    end
    if (strcmp(yname,'error_pAbs') || strcmp(yname,'error_p') || strcmp(yname,'surfaceError') || strcmp(yname,'energyError')) && (options.yScale == 1)
        yLabel = [yLabel, ' [\%]'];
    end
            
            
    xlabel(xLabel,'interpreter',intrprtrX)
    ylabel(yLabel,'interpreter',intrprtrY)
    h = gca;
    noGraphs = numel(h.Children);
    col = LaTeXcolorMap(noGraphs);
    for i = 1:noGraphs
        h.Children(i).Color = col(i,:);
    end
    savefig([subFolderName '/_' model '_' yname 'VS' xname])
end



