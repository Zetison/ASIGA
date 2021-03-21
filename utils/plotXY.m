function h = plotXY(x,y,axisType,lineStyle,color,legendEntry)
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