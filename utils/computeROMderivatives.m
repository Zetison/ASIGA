function task = computeROMderivatives(task,shiftROM)
% Construct the subspace V (scaled version of U = Pinv*V*omega^j)
if nargin < 2
    shiftROM = 0;
end
if task.misc.printLog
    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Computing derivatives for ROM ... ')
end

omega = task.misc.omega;
noRHSs = task.noRHSs;
dAdomega = cell(1,4);
dAdomega{1} = task.A1 + 2*omega*task.A2 + 4*omega^3*task.A4;
dAdomega{2} = 2*task.A2 + 12*omega^2*task.A4;
dAdomega{3} = 24*omega*task.A4;
dAdomega{4} = 24*task.A4;
if shiftROM == 0
    task.U = zeros(size(task.FF));
else
    task.U = [task.U, zeros(size(task.FF))];
end
dA = decomposition(task.Pinv*task.A*task.Pinv,'lu');
for i = 1:noRHSs
    j = shiftROM + i-1;
    b = task.FF(:,i);
    for k = 1:min(j,numel(dAdomega))
        b = b - nchoosek(j,k)*dAdomega{k}*task.U(:,j-k+1);
    end
    task.U(:,shiftROM+i) = task.Pinv*(dA\(task.Pinv*b));
end
i = 1:noRHSs;
j = shiftROM + i-1;
% keyboard
if true
    task.V = [task.V, task.U(:,shiftROM+i)];
else
    rightV = sqrt(omega.^j./factorial(j));
    task.V = [task.V, task.U(:,shiftROM+i).*rightV];
end
% solidIdx = task.varCol{1}.Aindices{2};
% leftV = ones(task.varCol{1}.noRows_tot,1);
% leftV(solidIdx) = leftV(solidIdx)*omega^2*task.varCol{1}.rho;
% allDofsToRemove = task.varCol{1}.allDofsToRemove;
% leftV(allDofsToRemove) = [];
% 
% task.V(:,indices) = leftV.*task.V(:,indices);
% 
% task.leftV = spdiags([diag(task.leftV); leftV],0,size(task.V,1),size(task.V,1));
% task.leftVinv = spdiags(1./diag(task.leftV),0,size(task.V,1),size(task.V,1));
% task.rightV = spdiags([diag(task.rightV); rightV.'],0,size(task.V,2),size(task.V,2));

if task.misc.printLog
    fprintf('using %12f seconds.', toc)
end