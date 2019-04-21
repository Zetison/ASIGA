function y = interPade(x,a,p,q)
inter = (a(2:end)+a(1:end-1))/2;
noDofs = size(p{1},2);
y = zeros(noDofs,numel(x),class(x));
for i = 1:numel(a)
    if i == 1
        indices = x <= inter(1);
    elseif i == numel(a)
        indices = x > inter(end);
    else
        indices = and(inter(i-1) < x, x <= inter(i));
    end
%     for j = 1:noDofs
%         y(j,indices) = polyval(p{i}(:,j),x(indices))./polyval(q{i}(:,j),x(indices));
%     end
    X = x(indices);
    P = size(p{i},1);
    Q = size(q{i},1);
    PQ = max(P,Q);
    XX = zeros(PQ,numel(X));
    for j = 1:PQ
        XX(j,:) = (X-a(i)).^(j-1); % p(1)*X^N + p(2)*X^(N-1) + ... + p(N)*X + p(N+1)
    end
    y(:,indices) = (p{i}.'*XX(1:P,:))./(q{i}.'*XX(1:Q,:));
end