function x_extremas = findExtremas(f, a, b, N)
newEpsilon = eps;
x = linspace(a, b, N);

tic
F = f(x);
toc

candidates = or(F(2:end-1) > max(F(1:end-2),F(3:end)), F(2:end-1) < min(F(1:end-2),F(3:end)));
candidates = find([0 candidates 0]);

x = [x(candidates-1); x(candidates+1)];
x_extremas = zeros(size(x));
options = optimset('TolX',newEpsilon, 'TolFun', newEpsilon, 'MaxIter', 1000,'MaxFunEvals',1000);%,'Display','iter');
parfor i = 1:size(x,2)
    tic
    I = x(:,i);
    x0 = mean(I);
    x_min = fminsearchbnd(f, x0, I(1), I(2), options);
    x_max = fminsearchbnd(@(x)-f(x), x0, I(1), I(2), options);
    x_res = zeros(2,1);
    if min(abs(x_min-I)) > 10*newEpsilon
        x_res(1) = x_min;
    end
    if min(abs(x_max-I)) > 10*newEpsilon
        x_res(2) = x_max;
    end
    x_extremas(:,i) = x_res;
    fprintf('Completed %d out of %d. Ellapsed time = %f\n', i, size(x,2), toc);
end
x_extremas = [x_extremas(1,:) x_extremas(2,:)];
x_extremas(x_extremas == 0) = [];