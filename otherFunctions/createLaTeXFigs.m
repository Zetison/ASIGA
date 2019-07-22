function createLaTeXFigs(newOptions)

options = struct('axisType',        'plot', ... 
                 'lineStyle',       '*-', ... 
                 'xLim',            [-inf,inf], ...
                 'yLim',            [-inf,inf], ...
                 'xLabel',          'x', ...
                 'yLabel',          'y', ...
                 'xScale',          1, ...
                 'enlargeXlimits',  false, ...
                 'enlargeYlimits',  false, ...
                 'yScale',          1, ... 
                 'colors',          [], ...
                 'legendEntries',   {{}}, ...
                 'legendLocation',  'northwest', ...
                 'xticks',          [], ...
                 'yticks',          [], ...
                 'xticklabels',     {{}}, ...
                 'xtickFormat',     '%g', ...
                 'addSlopes',       false, ...
                 'addSlopesExep',   [], ...
                 'fileNames',       {{}});
             
if nargin > 0
    newOptionFields = fieldnames(newOptions);
    for j = 1:numel(newOptionFields)
        options.(newOptionFields{j}) = newOptions.(newOptionFields{j});
    end
end
set(0,'defaultTextInterpreter','latex'); %trying to set the default
fileNames = options.fileNames;
legendEntries = options.legendEntries;
axisType = options.axisType;
M = length(fileNames);
if isempty(options.colors)
    col = LaTeXcolorMap(M);
else
    col = options.colors;
end
xMin = Inf;
xMax = -Inf;
yMin = Inf;
yMax = -Inf;

for i = 1:M
    X = readLaTeXFormat(fileNames{i});
    if iscell(options.lineStyle)
        lineStyle = options.lineStyle{i};
    else
        lineStyle = options.lineStyle;
    end
    h(i) = plotXY(X(:,1),X(:,2),axisType,lineStyle,col(i,:),'');
    if min(X(:,1)) < xMin
        xMin = min(X(:,1));
    end
    if max(X(:,1)) > xMax
        xMax = max(X(:,1));
    end
    if min(X(:,2)) < yMin
        yMin = min(X(:,2));
    end
    if max(X(:,2)) > yMax
        yMax = max(X(:,2));
    end
    hold on
    if options.addSlopes && ~any(i == options.addSlopesExep)
        x = log10(X(end-1,:));
        y = log10(X(end,:));
        x2 = x + 0.1*(y-x);
        y2 = x + 0.9*(y-x);
        if false
            plot(10.^[x2(1),x2(1),y2(1)],10.^[x2(2),y2(2),y2(2)],'black')
            text(10.^(x2(1)-0.02),10.^((x2(2)+y2(2))/2),num2str((y2(2)-x2(2))/(y2(1)-x2(1)),3),'HorizontalAlignment', 'right');
        else
            plot(10.^[x2(1),y2(1),y2(1)],10.^[x2(2),x2(2),y2(2)],'black')
            text(10.^(y2(1)+0.02),10.^((x2(2)+y2(2))/2),num2str((y2(2)-x2(2))/(y2(1)-x2(1)),3),'HorizontalAlignment', 'left');
        end
    end
end
if ~isempty(options.xticks)
    if strcmp(axisType,'polar')
        thetaticks(options.xticks)
    else
        xticks(options.xticks)
    end
end
if ~isempty(options.yticks)
    if strcmp(axisType,'polar')
        rticks(options.yticks)
    else
        yticks(options.yticks)
    end
end
if ~isempty(options.xticklabels)
    xticklabels(options.xticklabels)
end
if options.enlargeXlimits
    xLim = [xMin,xMax];
    switch axisType
        case {'loglog', 'semilogx'}
            xLim = log10(xLim);
            diff = xLim(2)-xLim(1);
            xlim(10.^(xLim + 0.1*[-1,1]*diff))
        otherwise
            diff = xLim(2)-xLim(1);
            xlim(xLim + 0.1*[-1,1]*diff)
    end
else
    if strcmp(axisType,'polar')
        if isinf(options.xLim(1)) && isinf(options.xLim(2))
            thetalim([xMin,xMax])
        else
            thetalim(options.xLim)
        end
    else
        if isinf(options.xLim(1)) && isinf(options.xLim(2))
            xlim([xMin,xMax])
        else
            xlim(options.xLim)
        end
    end
end
if options.enlargeYlimits
    yLim = [yMin,yMax];
    switch axisType
        case {'loglog', 'semilogy'}
            yLim = log10(yLim);
            diff = yLim(2)-yLim(1);
            ylim(10.^(yLim + 0.1*[-1,1]*diff))
        case 'polar'
            yLim = log10(yLim);
            diff = yLim(2)-yLim(1);
            rlim(10.^(yLim + 0.1*[-1,1]*diff))
        otherwise
            diff = yLim(2)-yLim(1);
            ylim(yLim + 0.1*[-1,1]*diff)
    end
else
    if strcmp(axisType,'polar')
        if isinf(options.yLim(1)) && isinf(options.yLim(2))
            rlim([yMin,yMax])
        else
            rlim(options.yLim)
        end
    else
        if isinf(options.yLim(1)) && isinf(options.yLim(2))
            ylim([yMin,yMax])
        else
            ylim(options.yLim)
        end
    end
end

if strcmp(axisType,'polar')
    thetatickformat(options.xtickFormat)
else
    xtickformat(options.xtickFormat)
    xlabel(options.xLabel,'interpreter','latex')
    ylabel(options.yLabel,'interpreter','latex')
end

if ~isempty(legendEntries)
    leg1 = legend(h(1:numel(legendEntries)),legendEntries);
    set(leg1,'Interpreter','latex','Location',options.legendLocation);
else
    legend off
end
