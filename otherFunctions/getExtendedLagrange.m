function V = getExtendedLagrange(X,n)
N = numel(X);
dp = zeros(N,N,N);
for i = 1:N
    dp(:,1,i) = lagrangePolynomials(X,i,N,X);
    dp(:,2:N,i) = lagrangePolynomialsNthDeriv(X,i,N,X,N-1);
end
diffr(1,1).coeffs = 1;
diffr(1,1).pows = [1,zeros(1,N-1)];
diffr(n,n).coeffs = 1;
diffr(n,n).pows = [1,zeros(1,N-1)];
for i = 1:n
    diffr(i,1).coeffs = 1;
    diffr(i,1).pows = [i,zeros(1,N-1)];
end
for k = 1:n
    for m = 2:n
        counter = 1;
        for i = 1:size(diffr(k,m-1).pows,1)
            pows = diffr(k,m-1).pows(i,:);
            coeffs = diffr(k,m-1).coeffs(i);
            indices = find(pows);
            for j = 1:numel(indices)
                idx = indices(j);
                exponent = pows(idx);
                powsTemp = pows;
                if idx < N
                    powsTemp(idx) = powsTemp(idx)-1;
                    powsTemp(idx+1) = powsTemp(idx+1)+1;
                    diffr(k,m).coeffs(counter) = coeffs*exponent;
                    diffr(k,m).pows(counter,:) = powsTemp;
                end
                counter = counter + 1;
            end
        end  
        [~, idxMap] = uniquetol(diffr(k,m).pows,eps,'ByRows',true, 'OutputAllIndices', true);    
        for i = 1:numel(idxMap)
            if numel(idxMap{i}) > 1
                diffr(k,m).coeffs(idxMap{i}(1)) = sum(diffr(k,m).coeffs(idxMap{i}));
                diffr(k,m).coeffs(idxMap{i}(2:end)) = [];
                diffr(k,m).pows(idxMap{i}(2:end),:) = [];
            end
        end
    end
end
V = zeros(n*N);
counter = 1;
for i = 1:n
    for j = 1:N
        temp = dp(j,:,:);
        temp = reshape(temp,N,N);
        temp2 = zeros(n,N);
        for m = 1:n
            for ii = 1:numel(diffr(m,i).coeffs)
                temp3 = 1;
                for jj = 1:size(diffr(m,i).pows,2)
                    temp3 = temp3.*temp(jj,:).^diffr(m,i).pows(ii,jj);
                end
                temp2(m,:) = temp2(m,:) + diffr(m,i).coeffs(ii)*temp3;
            end
        end
        V(counter,:) = reshape(temp2.',1,n*N);
        counter = counter + 1;
    end
end
        
end   
        

