function yq = interp1_new(x, y, xq)

[M,d,~] = size(y);
yq = zeros(M,d,length(xq));

y(:,:,end+1) = y(:,:,1);
for i = 1:M
    for j = 1:d
        y_temp = y(i,j,:);
        yq(i,j,:) = interp1(x, y_temp(:), xq, 'PCHIP');
    end
end