
a = get(gca,'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');
set(0,'DefaultLegendAutoUpdate','off')
for i = 1:numel(xdata)
    x = xdata{i};
    y = ydata{i};
    X = [x.',y.'];
    x1 = log10(X(end-1,:));
    y1 = log10(X(end,:));
    x2 = x1 + 0.1*(y1-x1);
    y2 = x1 + 0.9*(y1-x1);
    h = plot(10.^[x2(1),y2(1),y2(1)],10.^[x2(2),x2(2),y2(2)],'black');
    text(10.^(y2(1)+0.02),10.^((x2(2)+y2(2))/2),num2str((y2(2)-x2(2))/(y2(1)-x2(1)),4),'HorizontalAlignment', 'left');
end
set(0,'DefaultLegendAutoUpdate','on')