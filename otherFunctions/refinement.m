function [xNew, yNew, Fnew, refinementMap,totFuncEval] = refinement(x, y, f, F, delta)

xNew = zeros(1,2*length(x)-1);
yNew = zeros(1,2*length(y)-1);
Fnew = zeros(2*length(y)-1, 2*length(x)-1);
refinementMap = zeros(2*length(y)-1, 2*length(x)-1);
totFuncEval = 0;
for i = 1:length(x)-1
    xNew(2*i-1) = x(i);
    xNew(2*i) = (x(i)+x(i+1))/2;
end
xNew(end) = x(end);
for j = 1:length(y)-1
    yNew(2*j-1) = y(j);
    yNew(2*j) = (y(j)+y(j+1))/2;
end
yNew(end) = y(end);

%% Copy elements from F to Fnew
for i = 1:length(x)
    for j = 1:length(y)
        Fnew(2*j-1,2*i-1) = F(j,i);
    end
end

%% Calculate elements between columns
% Calculated elements cannot be directly inserted into Fnew due to the
% parfor loop. A tempMatrix is therefore made
tempMatrix = zeros(length(y), length(x)-1);
tempMatrix2 = zeros(length(y), length(x)-1);
xAverages = (x(1:end-1)+x(2:end))/2;
parfor i = 1:length(x)-1
    % Calculate the intermediate values
    temp_col = zeros(length(y),1);
    temp_col2 = zeros(length(y),1);
    xValue = xAverages(i);
    for j = 1:length(y)
        if abs(F(j,i) - F(j,i+1)) > delta
            temp_col(j) = f(xValue,y(j));
            temp_col2(j) = 1;
            totFuncEval = totFuncEval+1;
        else
            temp_col(j) = (F(j,i) + F(j,i+1))/2;
        end
    end
    tempMatrix(:,i) = temp_col;
    tempMatrix2(:,i) = temp_col2;
end

% Insert calculated elements into Fnew
for i = 1:length(x)-1
    for j = 1:length(y)
        Fnew(2*j-1,2*i) = tempMatrix(j,i);
    end
end
% Insert calculated elements into Fnew
for i = 1:length(x)-1
    for j = 1:length(y)
        refinementMap(2*j-1,2*i) = tempMatrix2(j,i);
    end
end

%% Calculate elements between rows
% Calculated elements cannot be directly inserted into Fnew due to the
% parfor loop. A tempMatrix is therefore made
tempMatrix = zeros(length(y)-1, length(x));
tempMatrix2 = zeros(length(y)-1, length(x));
yAverages = (y(1:end-1)+y(2:end))/2;
parfor i = 1:length(x)
    % Calculate the intermediate values
    temp_col = zeros(length(y)-1,1);
    temp_col2 = zeros(length(y)-1,1);
    xValue = x(i);
    for j = 1:length(yAverages)
        if abs(F(j,i) - F(j+1,i)) > delta
            temp_col(j) = f(xValue,yAverages(j));
            temp_col2(j) = 1;
        else
            temp_col(j) = (F(j,i) + F(j+1,i))/2;
        end
    end
    tempMatrix(:,i) = temp_col;
    tempMatrix2(:,i) = temp_col2;
end

% Insert calculated elements into Fnew
for i = 1:length(x)
    for j = 1:length(y)-1
        Fnew(2*j,2*i-1) = tempMatrix(j,i);
    end
end
% Insert calculated elements into Fnew
for i = 1:length(x)
    for j = 1:length(y)-1
        refinementMap(2*j,2*i-1) = tempMatrix2(j,i);
    end
end

%% Calculate elements on diagonals
% Calculated elements cannot be directly inserted into Fnew due to the
% parfor loop. A tempMatrix is therefore made
tempMatrix = zeros(length(y)-1, length(x)-1);
tempMatrix2 = zeros(length(y)-1, length(x)-1);
parfor i = 1:length(x)-1
    % Calculate the intermediate values
    temp_col = zeros(length(y)-1,1);
    temp_col2 = zeros(length(y)-1,1);
    xValue = xAverages(i);
    for j = 1:length(yAverages)
        if (abs(F(j,i) - F(j+1,i+1))/sqrt(2) > delta) || (abs(F(j,i+1) - F(j+1,i))/sqrt(2) > delta)
            temp_col(j) = f(xValue,yAverages(j));
            temp_col2(j) = 1;
        else
            temp_col(j) = (F(j,i) + F(j+1,i) + F(j,i+1) + F(j+1,i+1))/4;
        end
    end
    tempMatrix(:,i) = temp_col;
    tempMatrix2(:,i) = temp_col2;
end

% Insert calculated elements into Fnew
for i = 1:length(x)-1
    for j = 1:length(y)-1
        Fnew(2*j,2*i) = tempMatrix(j,i);
    end
end

% Insert calculated elements into Fnew
for i = 1:length(x)-1
    for j = 1:length(y)-1
        refinementMap(2*j,2*i) = tempMatrix2(j,i);
    end
end

