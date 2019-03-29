clear all
close all

% addpath ../integration

u = @(x) (1-x.^2).*exp(x);
dudx = @(x) (1-x.^2-2*x).*exp(x);
f = @(x) exp(x).*(x.^2+4*x+1);
% u = @(x) 1-x.^4;
% dudx = @(x) -4*x.^3;
% f = @(x) 12*x.^2;

% u = @(x) 1-x.^2;
% dudx = @(x) -2*x;
% f = @(x) 2*ones(size(x));

varCol.f = f;
varCol.analytic = u;
varCol.dudx = dudx;
options = {'operator','laplace',...
           'fieldDimension', 1,...
            'buildMassMatrix',0,...
           'applyBodyLoading', 1};
integrateRHSexact = 0;

N_arr = 3:33; % 33
E = zeros(size(N_arr));
cnd = zeros(size(N_arr));
E_IGA = zeros(size(N_arr));
cnd_IGA = zeros(size(N_arr));
E_IGA2 = zeros(size(N_arr));
cnd_IGA2 = zeros(size(N_arr));
noDofs_arr = zeros(size(N_arr));
for N_i = 1:numel(N_arr)
    N = N_arr(N_i);
    [tXi,rho] = GLLpoints(N);
    indices = 2:N-1;
    if integrateRHSexact
        nxi = N+10;
        [Q,W] = GLLpoints(nxi);
        Bxi = zeros(N,nxi);

        for i = 1:N
            Bxi(i,:) = lagrangePolynomials(Q,i,N,tXi);
        end
        b = Bxi*(f(Q).*W).';
    else
        b = (rho.*f(tXi)).';
    end
    D = derivativeMatrix(tXi,N);
    A = zeros(N-2);
    for i = indices
        for j = indices
            A(i-1,j-1) = sum(rho.*D(i,:).*D(j,:));
        end
    end
    U = A\b(indices);
    E(N_i) = computeH1error([0;U;0],tXi,N,varCol);
    cnd(N_i) = cond(A);
    
    nurbs = getRodData(-1,2);
    nurbs = elevateNURBSdegree(nurbs,N-2);
    varCol.dimension = 1;
    varCol = convertNURBS(nurbs, varCol);
    varCol = generateIGA1DMesh(varCol);
    [A, ~, F] = buildGlobalMatrices1D(varCol, options);
    dofsToRemove = [1,N];
    A(dofsToRemove,:) = [];
    A(:,dofsToRemove) = [];
    F(dofsToRemove) = [];
    U_IGA = A\F;
    E_IGA(N_i) = calcH1seminorm1D(varCol, [0;U_IGA;0]);
    cnd_IGA(N_i) = cond(full(A));
    
    
    nurbs = getRodData(-1,2);
    nurbs = elevateNURBSdegree(nurbs,N-2);
    nurbs = insertKnotsInNURBS(nurbs,linspace2(0,1,4));
    varCol.dimension = 1;
    varCol = convertNURBS(nurbs, varCol);
    varCol = generateIGA1DMesh(varCol);
    [A, ~, F] = buildGlobalMatrices1D(varCol, options);
    noDofs = varCol.patches{1}.noDofs;
    dofsToRemove = [1,noDofs];
    A(dofsToRemove,:) = [];
    A(:,dofsToRemove) = [];
    F(dofsToRemove) = [];
    U_IGA2 = A\F;
    E_IGA2(N_i) = calcH1seminorm1D(varCol, [0;U_IGA2;0]);
    cnd_IGA2(N_i) = cond(full(A));
    noDofs_arr(N_i) = noDofs-2;
end
figure(1)
noDofsSEM = N_arr-2;
semilogy(noDofsSEM,E,noDofsSEM,E_IGA,noDofs_arr,E_IGA2)
legend('SEM','IGA','IGA2')
printResultsToFile2('results/MA8502/poissonSEM', noDofsSEM.', E.')
printResultsToFile2('results/MA8502/poissonIGA1', noDofsSEM.', E_IGA.')
printResultsToFile2('results/MA8502/poissonIGA2', noDofs_arr.', E_IGA2.')
figure(2)
p = noDofsSEM.'+1;
p = linspace(p(1),p(end),1000).';
h = 2;
d = 1;
loglog(noDofsSEM,cnd, p-1, h^(-2)*p.^3,noDofsSEM,cnd_IGA, p-1, condBehaviorIGA(p,h,1),noDofs_arr,cnd_IGA2, p-1+4, condBehaviorIGA(p,h/5,1))
legend('SEM','SEM behavior','IGA','IGA behavior','IGA2','IGA2 behavior')
printResultsToFile2('results/MA8502/poissonSEM_cond', noDofsSEM.', cnd.')
printResultsToFile2('results/MA8502/poissonSEM_condBehavior', p-1, h^(-2)*p.^3)
printResultsToFile2('results/MA8502/poissonIGA1_cond', noDofsSEM.', cnd_IGA.')
printResultsToFile2('results/MA8502/poissonIGA1_condBehavior', p-1, condBehaviorIGA(p,h,1))
printResultsToFile2('results/MA8502/poissonIGA2_cond', noDofs_arr.', cnd_IGA2.')
printResultsToFile2('results/MA8502/poissonIGA2_condBehavior', p-1+4, condBehaviorIGA(p,h/5,1))


if 0
%     figure(42)
%     x = linspace(-1,1,1000).';
%     u_N = zeros(size(x));
%     for i = 2:N-1
%         u_N = u_N + U(i-1)*lagrangePolynomials(x,i,N,tXi);
%     end
%     plot(x,u(x))
%     hold on
%     plot(x,u_N)
%     legend('u','u_N')
    
    x = linspace(-1,1,1000).';
    close all
    N = 5;
%     tXi = GLLpoints(N);
    tXi = linspace(-1,1,N);
    figure(43)
    for i = 2
        p = lagrangePolynomials(x,i,N,tXi);
        plot(x,p,'DisplayName','dpdx0')
        dp = lagrangePolynomialsNthDeriv(x,i,N,tXi,3);
%         dp = lagrangePolynomialsNthDerivOld(x,i,N,tXi,2);
        hold on
        plot(x,dp(:,1),'DisplayName','dpdx1')
        hold on
        plot(x,dp(:,2),'DisplayName','dpdx2')
        hold on
        plot(x,dp(:,3),'DisplayName','dpdx3')
    end
    legend show
%     for i = 1:N
%         plot([tXi(i),tXi(i)],ylim,'--')
%     end
    plot(xlim,[0,0],'--','color','black')
%     
%     close all
%     N = 6;
%     tXi = linspace(-1,1,N);
%     tXi = GLLpoints(N);
%     for i = 1:N
%         for n = 1:20
%             p = (lagrangePolynomials(x,i,N,tXi)).^n;
%             hold on
%             plot(x,p,'DisplayName',['n = ' num2str(n) ', i = ' num2str(i)])
%         end
%     end
%     figure(44)
%     for i = 1:N
%         dpdxi = lagrangePolynomialsDeriv(x,i,N,tXi);
%         hold on
%         plot(x,dpdxi)
%     end
%     plot([x(1),x(end)],[0,0], 'Color', [17 17 17]/255)
%     grid on
end

return
error = 0;
N_arr = 3:200;
for N = N_arr
    tXi = lglnodes(N-1); 
    D = derivativeMatrix(tXi,N);
    d = diag(D);
    error = error + norm([d(1)+(N-1)*N/4; d(2:N-1); d(N)-(N-1)*N/4])/norm(d);
%     (N-1)*N/4
end
error

function K = condBehaviorIGA(p,h,d)

K = zeros(size(p));
indices = h < exp(-d*p/2);
K(indices) = h^(-2)*p(indices);
indices = and(exp(-d*p/2) < h, h < 1./p);
K(indices) = p(indices).*exp(d*p(indices));
indices = h > 1./p;
K(indices) = (exp(1)/4)^(d/h)*p(indices).^(-d/2).*h^(-d/2-1).*4.^(d*p(indices));

end




function relH1sError = computeH1error(U,tXi,N,varCol)

nxi = N+10;
[Q,W] = GLLpoints(nxi);
dBxi = zeros(N,nxi);

for i = 1:N
    dBxi(i,:) = lagrangePolynomialsDeriv(Q,i,N,tXi);
end

du_h = U.'*dBxi;

H1sError = 0;
normalizationH1s = 0;
dp = varCol.dudx(Q);
dp_e = dp-du_h;

dp2 = sum(abs(dp).^2,2);
dp_e2 = sum(abs(dp_e).^2,2);

H1sError = H1sError + sum(dp_e2.*W);
normalizationH1s = normalizationH1s + sum(dp2.*W);

relH1sError = 100*sqrt(H1sError/normalizationH1s);

end

