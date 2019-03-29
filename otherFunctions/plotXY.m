function h = plotXY(x,y,axisType,lineStyle,color,legendEntry)

switch axisType
    case 'plot'
        h = plot(x,y,lineStyle,'color',color,'DisplayName',legendEntry);
    case 'semilogx'
        h = semilogx(x,y,lineStyle,'color',color,'DisplayName',legendEntry);
    case 'semilogy'
        h = semilogy(x,y,lineStyle,'color',color,'DisplayName',legendEntry);
    case 'loglog'
        h = loglog(x,y,lineStyle,'color',color,'DisplayName',legendEntry);
end
legend('off');
legend('show');