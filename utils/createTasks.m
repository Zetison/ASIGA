function [tasks, counter, task] = createTasks(tasks, counter, task, i, loopParameters, loopParametersArr, connectedParameters, childrenParameters)
noParms = numel(loopParametersArr);
if i > noParms
    tasks(counter,1).task = task;
    counter = counter + 1;
else
    fieldName = loopParameters{noParms-i+1};
    temp = loopParametersArr{noParms-i+1};
    for j = 1:length(temp)
        if iscell(temp)
            eval(['task.' fieldName ' = temp{j};'])
        else
            eval(['task.' fieldName ' = temp(j);'])
        end
        for ii = 1:numel(connectedParameters)
            if isempty(connectedParameters{ii})
                continue
            end
            if strcmp(fieldName, connectedParameters{ii}{1})
                for jj = 2:numel(connectedParameters{ii})
                    temp2 = childrenParameters{ii}{jj};
                    if iscell(temp2)
                        eval(['task.' connectedParameters{ii}{jj} ' = temp2{j};'])
                    else
                        eval(['task.' connectedParameters{ii}{jj} ' = temp2(j);'])
                    end
                end
            end
        end
            
        [tasks, counter, task] = createTasks(tasks, counter, task, i+1, loopParameters, loopParametersArr, connectedParameters, childrenParameters);        
    end
end

end


