function mypostcallback(obj,evd,hobj1)
% keyboard
xLim = evd.Axes.XLim;
yLim = evd.Axes.YLim;
prev_xlim = hobj1.data.xLim;

x = hobj1.data.x;
y = hobj1.data.y;
f = hobj1.data.f;
if isrow(x)
    x = x';
end
if isrow(y)
    y = y';
end
originalXlim = hobj1.data.originalXlim;

specialValues = hobj1.data.specialXvalues;
noNewPoints = hobj1.data.noNewPoints;

if xLim(1) > prev_xlim(1) && xLim(2) < prev_xlim(2)
    indices = find((xLim(1) <= x).*(x <= xLim(2)));
    if length(indices) < noNewPoints
        disp('Inwards zoom encountered: new data points have to be calculated!')
        
        specialValuesOverlap = specialValues(find((xLim(1) < specialValues).*(specialValues < xLim(2))));
        new_xvalues = linspace(xLim(1),xLim(2),noNewPoints)';
        new_xvalues = sort([specialValuesOverlap; new_xvalues]);
        
        x(indices) = [];
        x = [x(1:(indices(1)-1)); new_xvalues; x(indices(1):end)];
        
        y(indices) = [];
        y = [y(1:indices(1)-1); f(new_xvalues).'; y(indices(1):end)];
        
        hobj1.data.x = x;
        hobj1.data.y = y;
        disp(['Number of datapoints is now ' num2str(length(x))])
    else
        disp('Inwards zoom encountered: no new data points needed...')
    end
end

plot(x,y);
% ticklabelformat(gca,'x','%1.10f')

xlim(xLim)
if isequal(xLim, originalXlim)
    ylim([min(y) 1.2*max(y)])
else
    ylim(yLim)
end
% ylim([0 1])

