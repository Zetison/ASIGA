function y = interPade(x,a,p,q)
inter = (a(2:end)+a(1:end-1))/2;
noDofs = size(p{1},2);
y = zeros(noDofs,numel(x));
for i = 1:numel(a)
    if i == 1
        indices = x <= inter(1);
    elseif i == numel(a)
        indices = x > inter(end);
    else
        indices = and(inter(i-1) < x, x <= inter(i));
    end
    for j = 1:noDofs
        y(j,indices) = polyval(p{i}(:,j),x(indices))./polyval(q{i}(:,j),x(indices));
    end
end