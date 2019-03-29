function ders = Bspline_basisDers3(i, u, p, U, n)
% This routine compute the p+1 nonzero basis functions and corresponding
% derivatives at xi

% Input
%       i:      knot span index corresonding to xi
%       p:      the degree of the B-Spline/NURBS
%       xi:     the value for which we want to evaluate the Bspline
%       Xi:     an open knot vector of size n+p+1

% Output
%       N:      array of the p+1 B-spline functions evaluated at xi  
i = i - 1;
left = zeros(1,p+1);
right = zeros(1,p+1);
ndu = zeros(p+1,p+1);
a = zeros(2,p+1);
ndu(0+1,0+1) = 1;
for j = 1:p
    left(j+1) = u-U(i+1-j+1);
    right(j+1) = U(i+j+1)-u;
    saved = 0;
    for r = 0:j-1
        ndu(j+1,r+1) = right(r+1+1)+left(j-r+1);
        temp = ndu(r+1,j-1+1)/ndu(j+1,r+1);
        ndu(r+1,j+1) = saved+right(r+1+1)*temp;
        saved = left(j-r+1)*temp;
    end
    ndu(j+1,j+1) = saved;
end
ders = zeros(n+1,p+1);
for j = 0:p
    ders(0+1,j+1) = ndu(j+1,p+1);
end
for r = 0:p
    s1 = 0; 
    s2 = 1;
    a(0+1,0+1) = 1;
    for k = 1:n
        d = 0;
        rk = r-k;
        pk = p-k;
        if r >= k
            a(s2+1,0+1) = a(s1+1,0+1)/ndu(pk+1+1,rk+1);
            d = a(s2+1,0+1)*ndu(rk+1,pk+1);
        end
        if rk >= -1
            j1 = 1;
        else
            j1 = -rk;
        end
        if r-1 <= pk
            j2 = k-1;
        else
            j2 = p-r;
        end
        for j = j1:j2
            a(s2+1,j+1) = (a(s1+1,j+1)-a(s1+1,j-1+1))/ndu(pk+1+1,rk+j+1);
            d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1);
        end
        if r <= pk
            a(s2+1,k+1) = -a(s1+1,k-1+1)/ndu(pk+1+1,r+1);
            d = d + a(s2+1,k+1)*ndu(r+1,pk+1);
        end
        ders(k+1,r+1) = d;
        j = s1;
        s1 = s2;
        s2 = j;
    end
end
r = p;
for k = 1:n
    for j = 0:p
        ders(k+1,j+1) = ders(k+1,j+1)*r;
    end
    r = r*(p-k);
end
        
% 
% varargout = cell(nargout,1);
% computeDers = nargout > 1;
% noxi = numel(xi);
% N = zeros(noxi,p+1);
% N(:,1) = 1;
% saved = ones(noxi,1);
% N_tilde = cell(p,1);
% N_tilde{end} = 1;
% for j = 2:p+1
%     % For k = 1 there is no dependence on N(k-1) of the previous run.
%     for k = 1:j
%         % Compute N_{i-j+k,j-1} according to the Cox-deBoor formula
%         temp = zeros(noxi,1);
%         if k ~= j
%             temp = (Xi(i+k)-xi)/(Xi(i+k)-Xi(i-j+k+1)).*N(:,k);
%         end
%         if k ~= 1
%             temp = temp + (xi-Xi(i-j+k))/(Xi(i+k-1)-Xi(i-j+k)).*saved;
%         end
%         saved = N(:,k);
%         N(:,k) = temp;
%     end
%     if j < p+1
%         N_tilde{p-j+1} = N(:,1:end-(p-j+1));
%     end
% %     if j == p
% %         N_tilde = N(:,1:end-1);
% %     elseif j == p-1
% %         N_tilde2 = N(:,1:end-2);
% %     elseif j == p-2
% %         N_tilde3 = N(:,1:end-3);
% %     end
% end
% % if p == 3
% %     N_tilde3 = 1;
% % end
% % 
% % if p == 2
% %     N_tilde2 = 1;
% % end
% % 
% % if p == 1
% %     N_tilde = 1;
% % end
% varargout{1} = N;
% if computeDers
%     alpha = cell(p+1,1);
%     alpha{1} = 1;
%     for k = 1:p
%         alpha{k+1} = zeros(1,k+1);
%         alpha{k+1}(1) = alpha{k}(1)./(Xi(i+1)-Xi(i-p+k));
%         alpha{k+1}(2:end-1) = (alpha{k}(2:end)-alpha{k}(1:end-1))./(Xi(i+2:i+p-k)-Xi(i-p+k+1:i-1));
%         alpha{k+1}(end) = -alpha{k}(end)./(Xi(i+p+1-k)-Xi(i));
%     end
%     
%     %% Calculate first derivatives
% %     dNdxi = zeros(noxi,p+1);
% % %     a = 1./(Xi(i+1:i+p)-Xi(i-p+1:i));
% % %     a = repmat(a,noxi,1);
% %     dNdxi(:,1:p)   = -N_tilde{1}.*a{1};
% %     dNdxi(:,2:p+1) = dNdxi(:,2:p+1) + N_tilde{1}.*a{1};
% %     varargout{2} = p*dNdxi;
% 
%     for k = 1:p
%         varargout{k+1} = zeros(noxi,p+1);
%         for j = 0:k
%             varargout{k+1}(:,j+1:p-k+j+1) = varargout{k+1}(:,j+1:p-k+j+1) + N_tilde{k}.*alpha{k+1}(j+1);
%         end
%         varargout{k+1} = prod(p:(p-k+1))*varargout{k+1};
%     end
% %     %% Calculate second derivatives
% %     d2Ndxi2 = zeros(noxi,p+1);
% %     if p > 1
% % %         a2 = 1./(Xi(i+1:i+p-1)-Xi(i-p+2:i));
% % %         a2 = repmat(a2,noxi,1);
% %         d2Ndxi2(:,1:p-1) =                    N_tilde{2}.*a{1}(:,1:end-1).*a{2};
% %         d2Ndxi2(:,2:p)   = d2Ndxi2(:,2:p)   - N_tilde{2}.*a{2}.*(a{1}(:,1:end-1) + a{1}(:,2:end));
% %         d2Ndxi2(:,3:p+1) = d2Ndxi2(:,3:p+1) + N_tilde{2}.*a{1}(:,2:end).*a{2};
% %         varargout{3} = p*(p-1)*d2Ndxi2;
% %     end
% % 
% %     %% Calculate third derivatives
% %     if nargout == 4
% %         d3Ndxi3 = zeros(noxi,p+1);
% %         if p > 2
% % %             a3 = 1./(Xi(i+1:i+p-2)-Xi(i-p+3:i));
% % %             a3 = repmat(a3,noxi,1);
% % 
% %             d3Ndxi3(:,1:p-2) =                  - N_tilde{3}.*a{1}(:,1:end-2).*a{2}(:,1:end-1).*a{3}; % 4
% %             d3Ndxi3(:,2:p-1) = d3Ndxi3(:,2:p-1) + N_tilde{3}.*a{3}.*(a{1}(:,2:end-1).*a{2}(:,2:end) + (a{1}(:,1:end-2) + a{1}(:,2:end-1)).*a{2}(:,1:end-1)); % 3
% %             d3Ndxi3(:,3:p)   = d3Ndxi3(:,3:p)   - N_tilde{3}.*a{3}.*((a{1}(:,2:end-1) + a{1}(:,3:end)).*a{2}(:,2:end) +  a{1}(:,2:end-1).*a{2}(:,1:end-1)); % 2
% %             d3Ndxi3(:,4:p+1) = d3Ndxi3(:,4:p+1) + N_tilde{3}.*a{1}(:,3:end).*a{2}(:,2:end).*a{3}; % 1
% %             varargout{4} = p*(p-1)*(p-2)*d3Ndxi3;
% %         end
% %     end
% %     
%     a = cell(p,1);
%     for j = 1:p
%         a{j} = repmat(1./(Xi(i+1:i+p-j+1)-Xi(i-p+j:i)),noxi,1);
%     end
% %     for k = 1:p
% %         varargout{k+1} = zeros(noxi,p+1);
% %         for j = 1:k+1
% %             temp = zeros(1,k+1);
% % %             for s = 1:
% %             varargout{k+1}(:,j:p+j-k) = varargout{k+1}(:,j:p+j-k) + (-1)^(j+k+1)*N_tilde{k}.*repmat(a{k}.*temp,noxi,1);
% %         end
% %         varargout{k+1} = prod(p:(p-k+1))*varargout{k+1};
% %     end
%         
% end
% 
% 
% 
% 
% 
