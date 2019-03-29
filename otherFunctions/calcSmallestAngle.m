function [underResolution, newXs] = calcSmallestAngle(x,y,minTheta)
y = abs(y);
y = y/max(y);
maxX = max(x);
x = x/maxX;
newXs = [];
underResolution = false;
insertedPreviousIteration = false;
for i = 1:length(x)-2
    vec1 = [x(i)-x(i+1), y(i)-y(i+1)];
    vec2 = [x(i+2)-x(i+1), y(i+2)-y(i+1)];
    theta = acos(dot(vec1, vec2)/(norm(vec1)*norm(vec2)));
    if theta < minTheta
        underResolution = true;
        
        if ~insertedPreviousIteration
            newXs = [newXs, (x(i)+x(i+1))/2];
        end
        newXs = [newXs, (x(i+1)+x(i+2))/2];
        insertedPreviousIteration = true;
    else
        insertedPreviousIteration = false;
    end
end
    
newXs = newXs*maxX;