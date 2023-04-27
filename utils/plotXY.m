function h = plotXY(x,y,axisType,lineStyle,color,legendEntry,addSlopes)
if all(isnan(x)) || all(isnan(y))
    error('All values are NaN')
end
if isempty(get(gca, 'Children'))
    xlimOld = [inf,-inf];
    ylimOld = [inf,-inf];
else
    if strcmp(axisType,'polar')
        xlimOld = thetalim;
        ylimOld = rlim;
    else
        xlimOld = xlim;
        ylimOld = ylim;
    end
end
if isempty(legendEntry)
    legendEntry = '';
end
switch axisType
    case 'plot'
        h = plot(x,y,lineStyle,'color',color,'DisplayName',legendEntry);
    case 'semilogx'
        h = semilogx(x,y,lineStyle,'color',color,'DisplayName',legendEntry);
    case 'semilogy'
        h = semilogy(x,y,lineStyle,'color',color,'DisplayName',legendEntry);
    case 'loglog'
        h = loglog(x,y,lineStyle,'color',color,'DisplayName',legendEntry);
    case 'polar'
        h = polarplot(x*pi/180,y,lineStyle,'color',color,'DisplayName',legendEntry);
end
xLim = [min(min(x),xlimOld(1)),max(max(x),xlimOld(2))];
yLim = [min(min(y),ylimOld(1)),max(max(y),ylimOld(2))];
if strcmp(axisType,'polar')
    if xLim(1) ~= xLim(2)
        thetalim(xLim)
    end
    if yLim(1) ~= yLim(2)
        rlim(yLim)
    end
%     rlim('auto')
else
    if xLim(1) ~= xLim(2)
        xlim(xLim)
    end
%     if yLim(1) ~= yLim(2)
%         ylim(yLim)
%     end
    ylim('auto')
end
legend('off');
leg1 = legend('show');
set(leg1,'Interpreter','latex');


hold on
if addSlopes && numel(x) > 1
    X = [x,y];
    x = log10(X(end-1,:));
    y = log10(X(end,:));
    x2 = x + 0.1*(y-x);
    y2 = x + 0.9*(y-x);
    if false
        plot(10.^[x2(1),x2(1),y2(1)],10.^[x2(2),y2(2),y2(2)],'black')
        text(10.^(x2(1)-0.02),10.^((x2(2)+y2(2))/2),num2str((y2(2)-x2(2))/(y2(1)-x2(1)),'%.3f'),'HorizontalAlignment', 'right');
    else
        plot(10.^[x2(1),y2(1),y2(1)],10.^[x2(2),x2(2),y2(2)],'black')
        text(10.^(y2(1)+0.02),10.^((x2(2)+y2(2))/2),num2str((y2(2)-x2(2))/(y2(1)-x2(1)),'%.3f'),'HorizontalAlignment', 'left');
    end
end