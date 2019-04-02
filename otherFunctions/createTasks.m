function [tasks, counter, task] = createTasks(tasks, counter, task, i, loopParameters, loopParametersArr)
noParms = numel(loopParametersArr);
if i > noParms
    tasks(counter,1).task = task;
    counter = counter + 1;
else
    fieldName = loopParameters{noParms-i+1};
    temp = loopParametersArr{noParms-i+1};
    for j = 1:length(temp)
        if iscell(temp)
            task.(fieldName) = temp{j};
        else
            task.(fieldName) = temp(j);
        end
            
        [tasks, counter, task] = createTasks(tasks, counter, task, i+1, loopParameters, loopParametersArr);        
    end
end

end


