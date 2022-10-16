function y = interTaylor(x,a,dF,D)
inter = (a(2:end)+a(1:end-1))/2;
y = zeros(size(dF{1},1),numel(x));
for i = 1:numel(a)
    if i == 1
        indices = x <= inter(1);
    elseif i == numel(a)
        indices = x > inter(end);
    else
        indices = and(inter(i-1) < x, x <= inter(i));
    end
    if any(indices)
        for n = 0:D
            y(:,indices) = y(:,indices) + dF{i}(:,n+1)/factorial(n)*(x(indices)-a(i)).^n;
        end
    end
end