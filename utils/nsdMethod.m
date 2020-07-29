function I = nsdMethod(cps,f,g,dgdx,omega,noGLpts,plotPaths,useFreudQuad,task)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotPaths
    figure(42)
end
Eps = 1e4*eps;
maxItrs = 100;
I = 0;

js = findDirsNSD(cps,dgdx);
dirs = findDirsNSD2(js,cps,dgdx);

for i = 1:size(cps,2)
    xi = cps(1,i);
    r = cps(2,i);
    if useFreudQuad
        [Q1D, W1D] = gaussFreudQuad(noGLpts,r);
    else
        [Q1D, W1D] = gaussLaguerreQuad(noGLpts,r);
    end
%     if task
%         keyboard
%     end
    if i ~= 1
%         xi_eps = xi + 1e-1*abs(mean([xi,cps(1,i-1)]))*h_end;
%         h_xi = xi_eps;
        
        F = F_(xi,r,omega,f,g,dgdx,dirs(2*i-2),Q1D,W1D,maxItrs,Eps,plotPaths,useFreudQuad);        
        I = I - F;
%         if task
%             keyboard
%             b = cps(1,end);
%             if i == 2
%                 F_exact = exp(1i*omega*g(xi))*integral(@(p) -f(-sqrt(abs(xi)+1i*p)).*1i./(2*sqrt(abs(xi)+1i*p)).*exp(-omega*p),0,inf);
%             elseif i == 3
%                 F_exact = exp(1i*omega*g(xi))*integral(@(p) f(sqrt(abs(xi)+1i*p)).*1i./(2*sqrt(abs(xi)+1i*p)).*exp(-omega*p),0,inf);
%             end
%             Error = 100*abs(F - F_exact)/abs(F_exact)
%             integral(@(x) f(x).*exp(1i*omega*g(x)),0,1)
%         end
    end
    if i ~= size(cps,2)
%         xi_eps = 0.9*xi + 0.1*cps(1,i+1) + sign(imag(h_end))*sign(r-1.5)*1e-4*1i;
%         h_xi = xi_eps;
        
        F = F_(xi,r,omega,f,g,dgdx,dirs(2*i-1),Q1D,W1D,maxItrs,Eps,plotPaths,useFreudQuad);
        I = I + F;
%         if task
%             if i == 1
%                 F_exact = exp(1i*omega*g(xi))*integral(@(p) -f(-sqrt(abs(xi)+1i*p)).*1i./(2*sqrt(abs(xi)+1i*p)).*exp(-omega*p),0,inf);
%             elseif i == 2
%                 F_exact = exp(1i*omega*g(xi))*integral(@(p) f(sqrt(abs(xi)+1i*p)).*1i./(2*sqrt(abs(xi)+1i*p)).*exp(-omega*p),0,inf);
%             end
% %             integral(@(p) 1i/2*exp(-omega*p),0,inf)
%             1i/2/omega
%             Error = 100*abs(F - F_exact)/abs(F_exact)
%             keyboard
%         end
    end
%     h_end = h_xi(end)/abs(h_xi(end));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotPaths
    gstr = func2str(g);
    nppts = 1000;
    if task == 2
        xLims = [0,1.4];
        yLims = [0,0.8];
        p_arr = linspace(0,32/25,1000).';
        printResultsToFile2('../results/NSD/task2path1', real(sqrt(1i*p_arr)), imag(sqrt(1i*p_arr)));  
        p_arr = linspace(0,8*sqrt(41)/25,1000).';
        printResultsToFile2('../results/NSD/task2path2', real(sqrt(1i*p_arr+1)), imag(sqrt(1i*p_arr+1)));  
%         plot(real(sqrt(1i*p_arr)), imag(sqrt(1i*p_arr)),'--') 
%         plot(real(sqrt(1i*p_arr+1)), imag(sqrt(1i*p_arr+1)),'--')
    else
        xLims = xlim;
        xLims = xLims + [-1,1]*0.03*max(abs(xLims));
        yLims = ylim;
        yLims = yLims + [-1,1]*0.03*max(abs(yLims));
    end
%     xLims = [-1.2, 1.2];
%     yLims = [-1, 1];
    x = linspace(xLims(1),xLims(2),nppts);
    y = linspace(yLims(1),yLims(2),nppts);
    [X,Y] = meshgrid(x,y);
    Z = imag(1i*g(X+1i*Y));
    
    switch task
        case 2
            lvlList = linspace(min(min(abs(Z))),max(max(abs(Z))),20);
        case {1,4,6,12,15}
            lvlList = linspace(min(min(abs(Z))),max(max(abs(Z))),20);
        case 13
            lvlList = 10.^linspace(log10(1e3*min(min(abs(Z)))),log10(max(max(abs(Z)))),20);
        case {14,16,17}
            lvlList = 10.^linspace(log10(1e6*min(min(abs(Z)))),log10(max(max(abs(Z)))),20);
        otherwise    
            lvlList = 10.^linspace(log10(1e4*min(min(abs(Z)))),log10(max(max(abs(Z)))),20);
        
    end
    lvlList = [-fliplr(lvlList), lvlList];
    [~, h1] = contour(x,y,Z,'LevelList',lvlList);
%     set(h1,'parent',ax(2));
    uistack(h1, 'bottom');
    colorbar
    caxis([-0.3,1.6])
%     axis equal tight
    fstr = func2str(f);
%     title(['g(x) = ' gstr(5:end) ' f(x) = ' fstr(5:end) ', \omega = ' num2str(omega)])
    
    drawnow
    if false
        for i = 1:size(cps,2)-1
            for j = 0:1
                xi = cps(1,i+j);
                d = dirs(2*i-1+j);
    %             d = s*d;
                if j == 0
                    quivers2(xi,0,real(d),imag(d),xLims,yLims,'black')
                else
                    dx = xLims(2)-xLims(1);
                    dy = yLims(2)-yLims(1);
                    theta = atan2(imag(d),real(d));
                    if theta == pi/2
                        s = 0.1*dy;
                    else
                        s = 0.1/sqrt(1/dx^2+tan(theta)^2/dy^2)/abs(cos(theta));
                    end
                    quivers2(real(xi+s*d),imag(s*d),-real(d),-imag(d),xLims,yLims,'black')
                end
            end
        end
    end
    hold off
    fileName = '../results/NSD/paths';
    extraAxisOptions = {...
        'axis on top=true', ...
        'at={(0,0)}', ...
        'xlabel=$\Re z$', ...
        'ylabel=$\Im z$', ...
        'colorbar style={ylabel={$\Im(\imag g(z))$}}'};

    matlab2tikz([fileName '_1.tex'], 'height', '3.2094in', ...
        'extraAxisOptions', extraAxisOptions, 'relativeDataPath', '../results/NSD/')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function F = F_(x,r,omega,f,g,dgdx,dirs,Q1D,W1D,maxItrs,Eps,plotPaths,useFreudQuad)
if ~useFreudQuad
    alpha = -r/(1+r);
end
F = 0;
for gp = 1:size(W1D,1)
    w = W1D(gp);
    if useFreudQuad
        u = Q1D(gp);
        q = u^(r+1);
        fact = (r+1)*u^r*w;
    else
        q = Q1D(gp);
        fact = q^(-alpha)*w;
    end
    p = q/omega;
    if r == 0
        gdx1 = dgdx{1}(x);
        gdx2 = dgdx{2}(x);
        if gdx2 == 0
            h_x = x + 1i*p/gdx1;
        else
            h_x = x - (gdx1-sign(gdx1)*sqrt(gdx1^2+2*1i*gdx2*p))/gdx2;
        end
    else
        h_x = x + dirs*(factorial(r+1)*p/abs(dgdx{r+1}(x)))^(1/(r+1));
    end
    [h_x,i] = newtonsMethod(@(hx)g(hx)-g(x)-1i*p,dgdx{1},h_x,maxItrs,Eps);
%     i
    dh_x = 1i/dgdx{1}(h_x);
    F = F + f(h_x)*dh_x*fact;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotPaths
%         f(h_x)*dh_x
        plot(h_x,'*','color','red')
%         plot(h_x,'*','color','blue')
        hold on
        h_x
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
F = F*exp(1i*omega*g(x))/omega;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotPaths
    nppts = 30;
    h_xi_plot = zeros(1,nppts);
    if useFreudQuad
        p2 = Q1D([1,end]).^(r+1)/omega;
    else
        p2 = Q1D([1,end])/omega;
    end
    p_plot = linspace(p2(1),p2(2),nppts);
    h_x = x;
    for j = 1:length(p_plot)
        p = p_plot(j);
        if abs(h_x-x) < 0.5
            if r == 0
                gdx1 = dgdx{1}(x);
                gdx2 = dgdx{2}(x);
                if gdx2 == 0
                    h_x = x + 1i*p/gdx1;
                else
                    h_x = x - (gdx1-sign(gdx1)*sqrt(gdx1^2+2*1i*gdx2*p))/gdx2;
                end
            else
                h_x = x + dirs*(factorial(r+1)*p/abs(dgdx{r+1}(x)))^(1/(r+1));
            end
        end
        h_x = newtonsMethod(@(hx)g(hx)-g(x)-1i*p,dgdx{1},h_x,maxItrs,Eps);
        h_xi_plot(j) = h_x;      
    end
    plot(h_xi_plot,'red')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



