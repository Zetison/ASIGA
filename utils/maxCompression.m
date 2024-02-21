function f = maxCompression(v,J,stopForDebug)

f = zeros(size(v));
for qp = 1:size(J{1},1)
    J_qp_T = [J{1}(qp,:); J{2}(qp,:); J{3}(qp,:)]; % The transpose of the Jacobian evaluated at quadrature point number qp
    E = 0.5*(J_qp_T + J_qp_T.');
    [V,D] = eig(E);
    eigenvalues = diag(D);
    [~,i] = max(eigenvalues);
%     f(qp,:) = abs(eigenvalues(i))*V(:,i);
    f(qp,:) = abs(eigenvalues(i))*V(i,:);
%     f(qp,:) = [0,0,-1];
%     f(qp,:) = V*eigenvalues;
%     if stopForDebug && qp == size(J{1},1)
%         keyboard
%     end
end