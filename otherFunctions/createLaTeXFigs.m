function createLaTeXFigs(newOptions)

options = struct('axisType',        'plot', ... 
                 'lineStyle',       '*-', ... 
                 'xLim',            [-inf,inf], ...
                 'yLim',            [-inf,inf], ...
                 'xLabel',          'x', ...
                 'yLabel',          'y', ...
                 'xScale',          1, ...
                 'yScale',          1, ... 
                 'colors',          [], ...
                 'legendEntries',   {{}}, ...
                 'legendLocation',  'northwest', ...
                 'xticks',          [], ...
                 'xticklabels',     {{}}, ...
                 'xtickFormat',     '%g', ...
                 'addSlopes',       false, ...
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
for i = 1:M
    X = readLaTeXFormat(fileNames{i});
    if iscell(options.lineStyle)
        lineStyle = options.lineStyle{i};
    else
        lineStyle = options.lineStyle;
    end
    h(i) = plotXY(X(:,1),X(:,2),axisType,lineStyle,col(i,:),'');
    hold on
    if options.addSlopes
        x = log10(X(end-1,:));
        y = log10(X(end,:));
        x2 = x + 0.1*(y-x);
        y2 = x + 0.9*(y-x);
        if i == 3
            plot(10.^[x2(1),x2(1),y2(1)],10.^[x2(2),y2(2),y2(2)],'black')
            text(10.^(x2(1)-0.02),10.^((x2(2)+y2(2))/2),num2str((y2(2)-x2(2))/(y2(1)-x2(1)),3),'HorizontalAlignment', 'right');
        else
            plot(10.^[x2(1),y2(1),y2(1)],10.^[x2(2),x2(2),y2(2)],'black')
            text(10.^(y2(1)+0.02),10.^((x2(2)+y2(2))/2),num2str((y2(2)-x2(2))/(y2(1)-x2(1)),3),'HorizontalAlignment', 'left');
        end
    end
end
if ~isempty(options.xticks)
    xticks(options.xticks)
end
if ~isempty(options.xticklabels)
    xticklabels(options.xticklabels)
end
xlim(options.xLim)
ylim(options.yLim)
xtickformat(options.xtickFormat)

xlabel(options.xLabel,'interpreter','latex')
ylabel(options.yLabel,'interpreter','latex')
    
leg1 = legend(h(1:numel(legendEntries)),legendEntries);
set(leg1,'Interpreter','latex','Location',options.legendLocation);

