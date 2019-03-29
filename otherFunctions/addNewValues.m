function [x,y] = addNewValues(x,y,x_new,y_new)

counter = 1;

for i = 1:length(x)
    if counter > length(x_new)
        return
    end
    if x(i+1) > x_new(counter)
        x = [x(1:i) x_new(counter) x(i+1:end)];
        y = [y(1:i) y_new(counter) y(i+1:end)];
        counter = counter + 1;
    end
end
    