clear all
close all

% addpath ../integration

u = @(x) x.*(1-x).*exp(x);
dudx = @(x) -(x.^2+x-1).*exp(x);
f = @(x) exp(x).*x.*(x+3);
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

N_arr = 10; % 33
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
    jacobian = 1/2;
    for i = indices
        for j = indices
            A(i-1,j-1) = sum(rho.*D(i,:).*D(j,:));
        end
    end
    U = A\b(indices);
    E(N_i) = computeH1error([0;U;0],tXi,N,varCol);
end


if 0
    figure(42)
    x = linspace(0,1,1000).';
    u_N = zeros(size(x));
    for i = 2:N-1
        u_N = u_N + U(i-1)*lagrangePolynomials(x,i,N,tXi);
    end
    plot(x,u(x))
    hold on
    plot(x,u_N)
    legend('u','u_N')
    
    N = 50;
%     tXi = GLLpoints(N);
    tXi = linspace(-1,1,N);
    figure(43)
    for i = 1
        p = lagrangePolynomials(x,i,N,tXi);
        hold on
        plot(x,p)
    end
    for i = 1:N
        yLim = ylim;
        plot([tXi(i),tXi(i)],yLim,'--')
    end
%     figure(44)
%     for i = 1:N
%         dpdxi = lagrangePolynomialsDeriv(x,i,N,tXi);
%         hold on
%         plot(x,dpdxi)
%     end
%     plot([x(1),x(end)],[0,0], 'Color', [17 17 17]/255)
%     grid on
end