clear all
close all
addpath export_fig

if true
    x = linspace(0,1,101);
    y = x;
    [X, Y] = meshgrid(x,y);
    X = X';
    Y = Y';
else
    X = 0.375;
    Y = 0.625;
end
u = zeros(size(X));
N = 0;
sourcePts = [0.5, 0.5;
             0.75, 0.75;
             0.25, 0.5];
strength = [ 1;
             1;
            -2];
if true
    while true
        maxTerm = 0;
        n = N;
        for m = 0:N
            if false
                fact_x = (2*m+1)*pi;
                fact_y = (2*n+1)*pi;
                a_mn = 16/(fact_x^2+fact_y^2)/pi^2/(4*m*n+2*m+2*n+1);
            else
                fact_x = (m+1)*pi;
                fact_y = (n+1)*pi;
                a_mn = 0;
                for i = 1:size(sourcePts,1)
                    a_mn = a_mn + strength(i)*sin(fact_x*sourcePts(i,1))*sin(fact_y*sourcePts(i,2));
                end
                a_mn = 4/(fact_x^2+fact_y^2)*a_mn;        
            end
            uTerm = a_mn*sin(fact_x*X).*sin(fact_y*Y);
            uTermAbs = max(max(abs(uTerm)));
            if uTermAbs > maxTerm
                maxTerm = uTermAbs;
            end
            u = u + uTerm;
            n = n - 1;
        end
        N = N + 1;
        maxu = max(max(abs(u)));
        relDiff = maxTerm/maxu;
        if mod(N,100) == 0
            fprintf('reldiff = %20.15g, u = %0.15g\n', relDiff, maxTerm)
        end
        if relDiff < 1e-5
            break
        end
    end
else
    m = inf;
    n = 0;
    prevAllZero = false;
    while m > 2 && ~prevAllZero
        maxTerm = 0;
        if m < 3
            prevAllZero = true;
        else
            prevAllZero = false;
        end
        m = 0;
        relDiff = inf;
        maxutermPrev = inf;
        while relDiff > 1e-10
            if false
                fact_x = (2*m+1)*pi;
                fact_y = (2*n+1)*pi;
                a_mn = 16/(fact_x^2+fact_y^2)/pi^2/(4*m*n+2*m+2*n+1);
            else
                fact_x = (m+1)*pi;
                fact_y = (n+1)*pi;
                temp = 0;
                for i = 1:size(sourcePts,1)
                    temp = temp + strength(i)*sin(fact_x*sourcePts(i,1))*sin(fact_y*sourcePts(i,2));
                end
                a_mn = 4/(fact_x^2+fact_y^2)*temp;        
            end
            uTerm = a_mn*sin(fact_x*X).*sin(fact_y*Y);
            
            u = u + uTerm;
            maxu = max(max(abs(u)));
            maxuterm = max(max(abs(uTerm)));
            relDiff = (maxuterm + maxutermPrev)/maxu;
            m = m + 1;
            maxutermPrev = maxuterm;
        end
        n = n + 1;
        if mod(n,1) == 0
            fprintf('reldiff = %20.15g, u = %0.15g\n', relDiff, maxu)
        end
    end
end
surf(X,Y,u,'EdgeColor','none','LineStyle','none');
colormap(jet)
camlight
grid off
axis off
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
export_fig('exact', '-png', '-transparent', '-r200')
% N = 10;
% A = zeros(N+1,N+1);
% for m = 0:N
%     for n = 0:N
%         fact_x = (m+1)*pi;
%         fact_y = (n+1)*pi;
%         temp = 0;
%         for i = 1:size(sourcePts,1)
%             temp = temp + strength(i)*sin(fact_x*sourcePts(i,1))*sin(fact_y*sourcePts(i,2));
%         end
%         A(m+1,n+1) = 4/(fact_x^2+fact_y^2)*temp;  
%     end
% end
% surf(log10(abs(A)))
% N = 100;
% A = zeros(N+1,N+1);
% for m = 0:N
%     for n = 0:N
%         fact_x = (m+1)*pi;
%         fact_y = (n+1)*pi;
%         A(m+1,n+1) = 4/(fact_x^2+fact_y^2);  
%     end
% end
% surf(log10(abs(A)))
% 
% m = inf;
% n = 0;
% prevAllZero = false;
% while m > 2 && ~prevAllZero
%     maxTerm = 0;
%     if m < 3
%         prevAllZero = true;
%     else
%         prevAllZero = false;
%     end
%     m = 0;
%     relDiff = inf;
%     maxutermPrev = inf;
%     while relDiff > 1e-10
%         if false
%             fact_x = (2*m+1)*pi;
%             fact_y = (2*n+1)*pi;
%             a_mn = 16/(fact_x^2+fact_y^2)/pi^2/(4*m*n+2*m+2*n+1);
%         else
%             fact_x = (m+1)*pi;
%             fact_y = (n+1)*pi;
%             temp = 0;
%             for i = 1:size(sourcePts,1)
%                 temp = temp + strength(i)*sin(fact_x*sourcePts(i,1))*sin(fact_y*sourcePts(i,2));
%             end
%             a_mn = 4/(fact_x^2+fact_y^2)*temp;        
%         end
%         uTerm = a_mn*sin(fact_x*X).*sin(fact_y*Y);
% 
%         u = u + uTerm;
%         maxu = max(max(abs(u)));
%         maxuterm = max(max(abs(uTerm)));
%         relDiff = (maxuterm + maxutermPrev)/maxu;
%         m = m + 1;
%         maxutermPrev = maxuterm;
%     end
%     n = n + 1;
%     if mod(n,1) == 0
%         fprintf('reldiff = %20.15g, u = %0.15g\n', relDiff, maxu)
%     end
% end
% -0.073671353281503

% -0.17252(1)

% M = 0:100;
% N = 0:100;
% [MM, NN] = meshgrid(M,N);
% surf(MM,NN, log(16./(((2*MM+1)*pi).^2+((2*NN+1)*pi).^2)/pi^2./(4*MM.*NN+2*MM+2*NN+1)))


% 
% 
% function uTerm = addContribution(m,n,x,y)
% % f = @(x,y) -2*(x-x.^2+y-y.^2);
% % f = @(x,y) exp(x+y);
% 
% 
% fact_x = (m+1)*pi;
% fact_y = (n+1)*pi;
% % N = 64;
% % M = 64;
% % MN = M*N;
% % [W,Q] = gaussianQuadNURBS(M,N);
% % a_mn = 0;
% % J_2 = 0.25;
% % for gp = 1:MN
% %     xit = Q(gp,1);
% %     etat = Q(gp,2);
% %     wt = W(gp);
% %     
% %     xt = (xit+1)/2;
% %     yt = (etat+1)/2;
% %     
% %     a_mn = a_mn + f(xt,yt)*sin(fact_x*xt)*sin(fact_y*yt)*J_2*wt;
% % end
% % integrand = @(xt,yt) f(xt,yt).*sin(fact_x*xt).*sin(fact_y*yt);
% % a_mn = integral2(integrand,0,1,0,1);
% a_mn = ((-1)^(n+m) + (-1)^n + (-1)^m + 1)/(m*n + m + n + 1)/pi^2;
% % 4/(fact_x^2+fact_y^2)
% a_mn = -4/(fact_x^2+fact_y^2)*a_mn;
% uTerm = a_mn*sin(fact_x*x)*sin(fact_y*y);
% 
