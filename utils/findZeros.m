function x_roots = findZeros(f,I,N)

x_arr = linspace(I(1),I(2),N);
f_arr = f(x_arr);
indices = find(diff(sign(f_arr)));

x_roots = zeros(numel(indices),1);
for i = 1:numel(indices)
    idx = indices(i);
    x_roots(i) = fzero(f, x_arr(idx:idx+1));
end
% f_prev = f(x);
% noZerosFound = 0;
% counter = 1;
%     
% while x < x_max
%     x = x + dx;
%     f_val = f(x);
%     if sign(f_prev) ~= sign(f_val)
% %         x_roots(counter) = bisection(f, x-dx, x, floor(log(dx/10^(-13))/log(2)), x*10^(-13));
%         x_roots(counter) = fzero(f, [x-dx x]);
%         x = x_roots(counter) + dx;
%         counter = counter + 1;
%         noZerosFound = noZerosFound + 1;
%         f_val = f(x);
%     end
%     f_prev = f_val;
% end
% x_roots = x_roots(1:(counter-1));