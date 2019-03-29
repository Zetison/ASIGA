function ss = matlabToComsol(n)

a = linspace(0,1,n).'.^2;
ss = [];
for i = 1:numel(a)
    ss = [ss, num2str(a(i))];
    if i ~= numel(a)
        ss = [ss, ', '];
    end
end