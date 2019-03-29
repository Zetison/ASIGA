function I = nsdMethod_2D(cps,f,g,dgdx,dgdy,omega,noGLpts,plotPaths,useFreudQuad,task,nwpe0,noGpts_,Eps,maxItrs)

a = cps.a;
b = cps.b;
c = cps.c;
d = cps.d;
% noGptsF = @(nwpe) 40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotPaths
    figure(42)
end
% if task
%     gstr = func2str(g);
%     nppts = 1000;
%     xLims = [-1,1];
%     yLims = [-0.15, 0.15];
%     
%     x = linspace(xLims(1),xLims(2),nppts);
%     y = linspace(yLims(1),yLims(2),nppts);
%     [X,Y] = meshgrid(x,y);
%     Z = imag(1i*g(X+1i*Y,c*ones(size(X))));
% 
% %     if strcmp(gstr(5:end),'x')
%     maxZ = max(max(abs(Z)));   
%     lvlstart = -1.5;
%     lvlList = 10.^linspace(lvlstart,log10(maxZ),40);
%     lvlList = [-fliplr(lvlList), lvlList];
% %         contour(x,y,Z,'LevelList',lvlList)
% 
%     [~,h1] = contour(x,y,Z,40);
%     uistack(h1, 'bottom');
% %     end
%     fstr = func2str(f);
%     title(['g(x) = ' gstr(7:end) ' f(x) = ' fstr(7:end) ', \omega = ' num2str(omega)])
%     if task == 14
%         xlim(xLims)
%         ylim(yLims)
%     end
%     colorbar
%     drawnow
%     hold on
%     keyboard
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yx = cps.yx;
xy = cps.xy;
I = 0;
f_arr = 1500*omega/(2*pi);
U = @(xi,eta,p) acos(cos(xi)-1i*p/cos(eta)/2);
V = @(xi,eta,q) acos(cos(eta)-1i*q/cos(xi)/2);
exactPath = 1;
% if f_arr > 2.4e4 % 2.385e4
%     keyboard;
% end
% if omega == 1000
%     keyboard
% end
nwpe_eta = -Inf;
for eta = [c,d]
    xiE = cps.cp_c(1,:);
    nwpe_eta = max([nwpe_eta, max(abs(omega*(g(xiE(1:end-1),eta)-g(xiE(2:end),eta))/(2*pi)))]);
end
nwpe_xi = -Inf;
for xi = [a,b]
    etaE = cps.cp_a(1,:);
    nwpe_xi = max([nwpe_xi, abs(omega*(g(xi,etaE(1:end-1))-g(xi,etaE(2:end)))/(2*pi))]);
end
    
if yx && nwpe_eta > nwpe0
    eta = c;   
    r_c = cps.cp_a(2,1); 
    jsx = findDirsNSD(cps.cp_c,dgdx,dgdy,1,eta);
    dirx = findDirsNSD2(jsx,cps.cp_c,dgdx,dgdy,1,eta);
    for i = 1:size(cps.cp_c,2)
        xi = cps.cp_c(1,i);
        r = cps.cp_c(2,i);
        jsy = findDirsNSD(cps.cp_a,dgdx,dgdy,2,xi);
        if i ~= 1 && prevNSD  
            I = I - Fyx_(xi,eta,[r,r_c],f,g,dgdx,dgdy,dirx(2*i-2),jsy,omega,noGLpts,plotPaths,useFreudQuad,cps.cp_a,Eps,maxItrs,exactPath, @(p) U(xi,eta,p), @(xi,q) V(xi,eta,q));
        end
        if i ~= size(cps.cp_c,2)
            xiE = cps.cp_c(1,i:i+1);
            nwpe = abs(omega*(g(xiE(1),eta)-g(xiE(2),eta))/(2*pi));
            if nwpe > nwpe0
                I = I + Fyx_(xi,eta,[r,r_c],f,g,dgdx,dgdy,dirx(2*i-1),jsy,omega,noGLpts,plotPaths,useFreudQuad,cps.cp_a,Eps,maxItrs,exactPath, @(p) U(xi,eta,p), @(xi,q) V(xi,eta,q));
                prevNSD = true;
            else
                noGpts = noGpts_(nwpe);
                diry = findDirsNSD2(jsy,cps.cp_a,dgdx,dgdy,2,xi);
                I = I + Fyxg_(eta,0,f,g,dgdy,diry,xiE,omega,noGpts,noGLpts,useFreudQuad,Eps,maxItrs,exactPath, @(xi,q) V(xi,eta,q));
                prevNSD = false;
            end
        end
    end
elseif xy && nwpe_xi > nwpe0
    xi = a;
    r_a = cps.cp_c(2,1);
    jsy = findDirsNSD(cps.cp_a,dgdx,dgdy,2,xi);
    diry = findDirsNSD2(jsy,cps.cp_a,dgdx,dgdy,2,xi);
    for i = 1:size(cps.cp_a,2)
        eta = cps.cp_a(1,i);
        r = cps.cp_a(2,i);
        jsx = findDirsNSD(cps.cp_c,dgdx,dgdy,1,eta);
        if i ~= 1 && prevNSD
            I = I - Fxy_(xi,eta,[r_a,r],f,g,dgdx,dgdy,jsx,diry(2*i-2),omega,noGLpts,plotPaths,useFreudQuad,cps.cp_c,Eps,maxItrs,exactPath, @(eta,p) U(xi,eta,p), @(q) V(xi,eta,q));
        end
        if i ~= size(cps.cp_a,2)
            etaE = cps.cp_a(1,i:i+1);
            nwpe = abs(omega*(g(xi,etaE(1))-g(xi,etaE(2)))/(2*pi));
            if nwpe > nwpe0
                I = I + Fxy_(xi,eta,[r_a,r],f,g,dgdx,dgdy,jsx,diry(2*i-1),omega,noGLpts,plotPaths,useFreudQuad,cps.cp_c,Eps,maxItrs, exactPath, @(eta,p) U(xi,eta,p), @(q) V(xi,eta,q));
                prevNSD = true;
            else
                noGpts = noGpts_(nwpe);
                dirx = findDirsNSD2(jsx,cps.cp_c,dgdx,dgdy,1,eta);
                I = I + Fxyg_(xi,0,f,g,dgdx,dirx,etaE,omega,noGpts,noGLpts,useFreudQuad,Eps,maxItrs, exactPath, @(eta,p) U(xi,eta,p));
                prevNSD = false;
            end
        end
    end
else
    [W2D,Q2D] = gaussianQuadNURBS(noGpts_(nwpe_xi),noGpts_(nwpe_eta));
%     [W2D,Q2D] = gaussianQuadNURBS(64,64);
    x  = parent2ParametricSpace([a,b],  Q2D(:,1));
    y  = parent2ParametricSpace([c,d],  Q2D(:,2));
    J2 = 1/4*(b-a)*(d-c);
    fact = J2*W2D;
    I = sum(f(x,y).*exp(1i*omega*g(x,y)).*fact);
end
if false
    eta = c;
    R = 5;
    dmx = -[1,2,3];
    dmx = 2*dmx/norm(dmx);
    gamma = @(xi) dmx(1)*cos(xi)+dmx(2)*sin(xi);
    
    v_eta = @(xi,q) -1i*log(((g(xi,eta)+1i*q)/R + sqrt((g(xi,eta)+1i*q).^2/R^2-dmx(3)^2-gamma(xi).^2))./(gamma(xi)-1i*dmx(3)));
    dv_eta = @(xi,q) 1i./dgdy{1}(xi,v_eta(xi,q));
    F_Exact2 = integral2(@(xi,q)exp(1i*omega*g(xi,eta)).*f(xi,v_eta(xi,q)).*exp(-omega*q).*dv_eta(xi,q),xiE(1),xiE(2),0,inf)
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    xi = cps.cp_c(1,1);
    
    
    u_xi = @(p) -1i*log(((g(xi,eta)+1i*p)/R-dmx(3)*sin(eta) - sqrt(((g(xi,eta)+1i*p)/R-dmx(3)*sin(eta)).^2-(dmx(1)^2+dmx(2)^2)*cos(eta).^2))/(dmx(1)-1i*dmx(2))./cos(eta));
    du_xi = @(p) 1i./dgdx{1}(u_xi(p),eta);
    v_eta = @(xi,q) -1i*log(((g(xi,eta)+1i*q)/R + sqrt((g(xi,eta)+1i*q).^2/R^2-dmx(3)^2-gamma(xi).^2))./(gamma(xi)-1i*dmx(3)));
    dv_eta = @(xi,q) 1i./dgdy{1}(xi,v_eta(xi,q));
    F_Exact = exp(1i*omega*g(xi,eta))*integral2(@(p,q)f(u_xi(p),v_eta(u_xi(p),q)).*exp(-omega*(p+q)).*du_xi(p).*dv_eta(u_xi(p),q),0,inf,0,inf)

    xi = cps.cp_c(1,2);
    u_xi = @(p) -1i*log(((g(xi,eta)+1i*p)/R-dmx(3)*sin(eta) - sqrt(((g(xi,eta)+1i*p)/R-dmx(3)*sin(eta)).^2-(dmx(1)^2+dmx(2)^2)*cos(eta).^2))/(dmx(1)-1i*dmx(2))./cos(eta));
    du_xi = @(p) 1i./dgdx{1}(u_xi(p),eta);
    v_eta = @(xi,q) -1i*log(((g(xi,eta)+1i*q)/R + sqrt((g(xi,eta)+1i*q).^2/R^2-dmx(3)^2-gamma(xi).^2))./(gamma(xi)-1i*dmx(3)));
    dv_eta = @(xi,q) 1i./dgdy{1}(xi,v_eta(xi,q));
    F_Exact = exp(1i*omega*g(xi,eta))*integral2(@(p,q)f(u_xi(p),v_eta(u_xi(p),q)).*exp(-omega*(p+q)).*du_xi(p).*dv_eta(u_xi(p),q),0,inf,0,inf)

    xi = cps.cp_c(1,2);
    u_xi = @(p) -1i*log(((g(xi,eta)+1i*p)/R-dmx(3)*sin(eta) + sqrt(((g(xi,eta)+1i*p)/R-dmx(3)*sin(eta)).^2-(dmx(1)^2+dmx(2)^2)*cos(eta).^2))/(dmx(1)-1i*dmx(2))./cos(eta));
    du_xi = @(p) 1i./dgdx{1}(u_xi(p),eta);
    v_eta = @(xi,q) -1i*log(((g(xi,eta)+1i*q)/R + sqrt((g(xi,eta)+1i*q).^2/R^2-dmx(3)^2-gamma(xi).^2))./(gamma(xi)-1i*dmx(3)));
    dv_eta = @(xi,q) 1i./dgdy{1}(xi,v_eta(xi,q));
    F_Exact = exp(1i*omega*g(xi,eta))*integral2(@(p,q)f(u_xi(p),v_eta(u_xi(p),q)).*exp(-omega*(p+q)).*du_xi(p).*dv_eta(u_xi(p),q),0,inf,0,inf)

    xi = cps.cp_c(1,3);
    u_xi = @(p) -1i*log(((g(xi,eta)+1i*p)/R-dmx(3)*sin(eta) + sqrt(((g(xi,eta)+1i*p)/R-dmx(3)*sin(eta)).^2-(dmx(1)^2+dmx(2)^2)*cos(eta).^2))/(dmx(1)-1i*dmx(2))./cos(eta));
    du_xi = @(p) 1i./dgdx{1}(u_xi(p),eta);
    v_eta = @(xi,q) -1i*log(((g(xi,eta)+1i*q)/R + sqrt((g(xi,eta)+1i*q).^2/R^2-dmx(3)^2-gamma(xi).^2))./(gamma(xi)-1i*dmx(3)));
    dv_eta = @(xi,q) 1i./dgdy{1}(xi,v_eta(xi,q));
    F_Exact = exp(1i*omega*g(xi,eta))*integral2(@(p,q)f(u_xi(p),v_eta(u_xi(p),q)).*exp(-omega*(p+q)).*du_xi(p).*dv_eta(u_xi(p),q),0,inf,0,inf)
    
end
%         xi = a;
%         v_eta = @(q) 1-xi-sqrt((1-xi)^2+eta^2-2*(1-xi)*eta+2*1i*q);
%         u_xi = @(p,eta) (xi*(eta+2)+1i*p)./(2+eta);
%         dv_eta = @(q) -1i./sqrt((1-xi)^2+eta^2-2*(1-xi)*eta+2*1i*q);
%         du_xi = @(p,eta) 1i./(2+eta);
%         F_Exact = exp(1i*omega*g(xi,eta))*integral2(@(p,q)f(u_xi(p,v_eta(q)),v_eta(q)).*exp(-omega*(p+q)).*du_xi(p,v_eta(q)).*dv_eta(q),0,inf,0,inf);
%         abs(F_Exact-F)/F_Exact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotPaths
    keyboard
    hold on
    gstr = func2str(g);
    nppts = 1000;
    xLims = xlim;
    xLims = xLims + [-1,1]*0.1*max(abs(xLims));
    yLims = ylim;
    yLims = yLims + [-1,1]*0.1*max(abs(yLims));
    
    x = linspace(xLims(1),xLims(2),nppts);
    y = linspace(yLims(1),yLims(2),nppts);
    [X,Y] = meshgrid(x,y);
    Z = imag(1i*g(X+1i*Y,c));

%     if strcmp(gstr(5:end),'x')
%     maxZ = max(max(abs(Z)));   
%     lvlstart = -1.5;
%     lvlList = 10.^linspace(lvlstart,log10(maxZ),40);
%     lvlList = [-fliplr(lvlList), lvlList];
%         contour(x,y,Z,'LevelList',lvlList)

    contour(x,y,Z,40)
    fstr = func2str(f);
    title(['g(x) = ' gstr(7:end) ' f(x) = ' fstr(7:end) ', \omega = ' num2str(omega)])
    if task == 14
        xlim(xLims)
        ylim(yLims)
    end
    colorbar
    drawnow
    hold off
    figure(43)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dirs = findDirsNSD(cps.cp_d,dgdx,dgdy,1,d);
% for i = 1:size(cps.cp_d,1)
%     xi = cps.cp_d(i,1);
%     r = cps.cp_d(i,2);
%     if i ~= 1
%         F = Fyx_(xi,d,[r,0],f,g,dgdx,dgdy,dirx(2*i-2),diry,omega,noGLpts,plotPaths,useFreudQuad);     
%         I = I + F;
%     end
%     if i ~= size(cps.cp_d,1)
%         F = Fyx_(xi,d,[r,0],f,g,dgdx,dgdy,dirx(2*i-1),diry,omega,noGLpts,plotPaths,useFreudQuad);  
%         I = I - F;
%     end
% end
if yx && nwpe_eta > nwpe0
    eta = d;    
    r_d = cps.cp_b(2,1);
    jsx = findDirsNSD(cps.cp_d,dgdx,dgdy,1,eta);
    dirx = findDirsNSD2(jsx,cps.cp_d,dgdx,dgdy,1,eta);
    for i = 1:size(cps.cp_d,2)
        xi = cps.cp_d(1,i);
        r = cps.cp_d(2,i);
        jsy = findDirsNSD(cps.cp_b,dgdx,dgdy,2,xi);
        if i ~= 1 && prevNSD   
            I = I + Fyx_(xi,eta,[r,r_d],f,g,dgdx,dgdy,dirx(2*i-2),jsy,omega,noGLpts,plotPaths,useFreudQuad,cps.cp_b,Eps,maxItrs, exactPath, @(p) U(xi,eta,p), @(xi,q) V(xi,eta,q));
        end
        if i ~= size(cps.cp_d,2)
            xiE = cps.cp_d(1,i:i+1);
            nwpe = abs(omega*(g(xiE(1),eta)-g(xiE(2),eta))/(2*pi));
            if nwpe > nwpe0
                I = I - Fyx_(xi,eta,[r,r_d],f,g,dgdx,dgdy,dirx(2*i-1),jsy,omega,noGLpts,plotPaths,useFreudQuad,cps.cp_b,Eps,maxItrs, exactPath, @(p) U(xi,eta,p), @(xi,q) V(xi,eta,q));
                prevNSD = true;
            else
                noGpts = noGpts_(nwpe);
                diry = findDirsNSD2(jsy,cps.cp_b,dgdx,dgdy,2,xi);
                I = I - Fyxg_(eta,0,f,g,dgdy,diry,xiE,omega,noGpts,noGLpts,useFreudQuad,Eps,maxItrs, exactPath, @(xi,q) V(xi,eta,q));
                prevNSD = false;
            end
        end
    end
elseif xy && nwpe_xi > nwpe0
    xi = b;
    r_b = cps.cp_d(2,1);
    jsy = findDirsNSD(cps.cp_b,dgdx,dgdy,2,xi);
    diry = findDirsNSD2(jsy,cps.cp_b,dgdx,dgdy,2,xi);
    jsx = findDirsNSD(cps.cp_d,dgdx,dgdy,1,eta);
    for i = 1:size(cps.cp_b,2)
        eta = cps.cp_b(1,i);
        r = cps.cp_b(2,i);
        if i ~= 1 && prevNSD
            I = I + Fxy_(xi,eta,[r_b,r],f,g,dgdx,dgdy,jsx,diry(2*i-2),omega,noGLpts,plotPaths,useFreudQuad,cps.cp_d,Eps,maxItrs,exactPath, @(eta,p) U(xi,eta,p), @(q) V(xi,eta,q));
        end
        if i ~= size(cps.cp_b,2)
            etaE = cps.cp_b(1,i:i+1);
            nwpe = abs(omega*(g(xi,etaE(1))-g(xi,etaE(2)))/(2*pi));
            if nwpe > nwpe0
                I = I - Fxy_(xi,eta,[r_b,r],f,g,dgdx,dgdy,jsx,diry(2*i-1),omega,noGLpts,plotPaths,useFreudQuad,cps.cp_d,Eps,maxItrs,exactPath, @(eta,p) U(xi,eta,p), @(q) V(xi,eta,q));
                prevNSD = true;
            else
                noGpts = noGpts_(nwpe);
                dirx = findDirsNSD2(jsx,cps.cp_d,dgdx,dgdy,1,eta);
                I = I - Fxyg_(xi,0,f,g,dgdx,dirx,etaE,omega,noGpts,noGLpts,useFreudQuad,Eps,maxItrs,exactPath, @(eta,p) U(xi,eta,p));
                prevNSD = false;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotPaths
    figure(42)
    hold on
    gstr = func2str(g);
    nppts = 1000;
    xLims = xlim;
    xLims = xLims + [-1,1]*0.1*max(abs(xLims));
    yLims = ylim;
    yLims = yLims + [-1,1]*0.1*max(abs(yLims));
    
    x = linspace(xLims(1),xLims(2),nppts);
    y = linspace(yLims(1),yLims(2),nppts);
    [X,Y] = meshgrid(x,y);
    Z = imag(1i*g(X+1i*Y,d));

%     if strcmp(gstr(5:end),'x')
    maxZ = max(max(abs(Z)));   
    lvlstart = -1.5;
    lvlList = 10.^linspace(lvlstart,log10(maxZ),40);
    lvlList = [-fliplr(lvlList), lvlList];
%         contour(x,y,Z,'LevelList',lvlList)

    [~,h1] = contour(x,y,Z,40);
    uistack(h1, 'bottom');
%     end
    fstr = func2str(f);
    title(['g(x) = ' gstr(7:end) ' f(x) = ' fstr(7:end) ', \omega = ' num2str(omega)])
    if task == 14
        xlim(xLims)
        ylim(yLims)
    end
    colorbar
    drawnow
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = Fxy_(x,y,r,f,g,dgdx,dgdy,jsx,diry,omega,noGLpts,plotPaths,useFreudQuad,cp,Eps,maxItrs,exactPath,u_x_,v_y_)

if useFreudQuad
    [Q1D_p, W1D_p] = gaussFreudQuad(noGLpts,r(1));
    [Q1D_q, W1D_q] = gaussFreudQuad(noGLpts,r(2));
else
    [Q1D_p, W1D_p] = gaussLaguerreQuad(noGLpts,r(1));
    [Q1D_q, W1D_q] = gaussLaguerreQuad(noGLpts,r(2));
    alpha = -r./(1+r);
end

gdy0 = g(x,y);
if r(2) == 0
    gdy1 = dgdy{1}(x,y);
    gdy2 = dgdy{2}(x,y);
end
F = 0;
for j = 1:size(W1D_q,1)
    wt_q = W1D_q(j);
    if useFreudQuad
        qtt = Q1D_q(j);
        qt = qtt^(r(2)+1);
    else
        qt = Q1D_q(j);
    end
    q = qt/omega;
    if exactPath
        v_y = v_y_(q);
    else
        if r(2) == 0
            if gdy2 == 0
                v_y = y + 1i*q/gdy1;
            else
                v_y = y - (gdy1-sign(real(gdy1))*sqrt(gdy1^2+2*1i*gdy2*q))/gdy2;
            end
        else
            v_y = y + diry*(factorial(r(2)+1)*q/abs(dgdy{r(2)+1}(x,y)))^(1/(r(2)+1));
        end
        [v_y,noitr] = newtonsMethod(@(vy)g(x,vy)-gdy0-1i*q,@(vy)dgdy{1}(x,vy),v_y,maxItrs,Eps);
    %     noitr
    end
    dv_y = 1i/dgdy{1}(x,v_y);

    gdx0 = g(x,v_y);
    if r(1) == 0
        gdx1 = dgdx{1}(x,v_y);
        gdx2 = dgdx{2}(x,v_y);
    end
%     dirx = findDirsNSD2(jsx,cp,dgdx,dgdy,1,v_y);
    for i = 1:size(W1D_p,1)
        wt = W1D_p(i)*wt_q;
        if useFreudQuad
            ptt = Q1D_p(i);
            pt = ptt^(r(1)+1);
        else
            pt = Q1D_p(i);
        end
        
        p = pt/omega;
        if useFreudQuad
            fact = (r(1)+1)*(r(2)+1)*ptt^r(1)*qtt^r(2)*wt;
        else
            fact = pt^(-alpha(1))*qt^(-alpha(2))*wt;
        end
        if exactPath
            u_x = u_x_(v_y,p);
        else
            if r(1) == 0
                if gdx2 == 0
                    u_x = x + 1i*p/gdx1;
                else
                    u_x = x - (gdx1-sign(real(gdx1))*sqrt(gdx1^2+2*1i*gdx2*p))/gdx2;
                end
            else
                dirx = findDirsNSD2(jsx,cp,dgdx,dgdy,1,v_y);
                u_x = x + dirx(cp(1,:) == x)*(factorial(r(1)+1)*p/abs(dgdx{r(1)+1}(x,v_y)))^(1/(r(1)+1));
            end

            [u_x,noitr] = newtonsMethod(@(ux)g(ux,v_y)-gdx0-1i*p,@(ux)dgdx{1}(ux,v_y),u_x,maxItrs,Eps);
    %         noitr
        end
        du_x = 1i/dgdx{1}(u_x,v_y);
        
        F = F + f(u_x, v_y)*dv_y*du_x*fact;
    end
end
F = F*exp(1i*omega*g(x,y))/omega^2;

function F = Fxyg_(x,r,f,g,dgdx,dirx,etaE,omega,noGpts,noGLpts,useFreudQuad,Eps,maxItrs,exactPath,u_x_)

if useFreudQuad
    [Q1D_p, W1D_p] = gaussFreudQuad(noGLpts,r);
    [Q1D_q, W1D_q] = gaussQuad(noGpts);
else
    [Q1D_p, W1D_p] = gaussLaguerreQuad(noGLpts,r);
    [Q1D_q, W1D_q] = gaussQuad(noGpts);
    alpha = -r/(1+r);
end

F = 0;
J_2 = 0.5*(etaE(2)-etaE(1));
for j = 1:size(W1D_q,1)
    wt_q = W1D_q(j)*J_2;
    eta = parent2ParametricSpace(etaE,Q1D_q(j));

    gdx0 = g(x,eta);
    if r == 0
        gdx1 = dgdx{1}(x,eta);
        gdx2 = dgdx{2}(x,eta);
    end
    expFact = exp(1i*omega*gdx0);
    
    for i = 1:size(W1D_p,1)
        wt = W1D_p(i)*wt_q;
        if useFreudQuad
            ptt = Q1D_p(i);
            pt = ptt^(r+1);
        else
            pt = Q1D_p(i);
        end
        
        p = pt/omega;
        if useFreudQuad
            fact = (r+1)*ptt^r*wt;
        else
            fact = pt^(-alpha)*wt;
        end
        if exactPath
            u_x = u_x_(eta,p);
        else
            if r == 0
                if gdx2 == 0
                    u_x = x + 1i*p/gdx1;
                else
                    u_x = x - (gdx1-sign(real(gdx1))*sqrt(gdx1^2+2*1i*gdx2*p))/gdx2;
                end
            else
                u_x = x + dirx*(factorial(r+1)*p/abs(dgdx{r+1}(x,eta)))^(1/(r+1));
            end
            u_x = newtonsMethod(@(ux)g(ux,eta)-gdx0-1i*p,@(ux)dgdx{1}(ux,eta),u_x,maxItrs,Eps);
        end
        
        du_x = 1i/dgdx{1}(u_x,eta);
        
        F = F + expFact*f(u_x, eta)*du_x*fact;
    end
end
F = F/omega;


function F = Fyxg_(y,r,f,g,dgdy,diry,xiE,omega,noGpts,noGLpts,useFreudQuad,Eps,maxItrs,exactPath,v_y_)

if useFreudQuad 
    [Q1D_p, W1D_p] = gaussQuad(noGpts);
    [Q1D_q, W1D_q] = gaussFreudQuad(noGLpts,r);
else
    [Q1D_p, W1D_p] = gaussQuad(noGpts);
    [Q1D_q, W1D_q] = gaussLaguerreQuad(noGLpts,r);
    alpha = -r/(1+r);
end

F = 0;
J_2 = 0.5*(xiE(2)-xiE(1));
for i = 1:size(W1D_p,1)
    wt_p = W1D_p(i)*J_2;
    xi = parent2ParametricSpace(xiE,Q1D_p(i));
    
    gdy0 = g(xi,y);
    if r == 0
        gdy1 = dgdy{1}(xi,y);
        gdy2 = dgdy{2}(xi,y);
    end
    expFact = exp(1i*omega*gdy0);
    
    for j = 1:size(W1D_q,1)
        wt = W1D_q(j)*wt_p;
        if useFreudQuad
            qtt = Q1D_q(j);
            qt = qtt^(r+1);
        else
            qt = Q1D_q(j);
        end
        
        q = qt/omega;
        if useFreudQuad
            fact = (r+1)*qtt^r*wt;
        else
            fact = qt^(-alpha)*wt;
        end
        if exactPath
            v_y = v_y_(xi,q);
        else
            if r == 0
                if gdy2 == 0
                    v_y = y + 1i*q/gdy1;
                else
                    v_y = y - (gdy1-sign(real(gdy1))*sqrt(gdy1^2+2*1i*gdy2*q))/gdy2;
                end
            else
                v_y = y + diry*(factorial(r+1)*q/abs(dgdy{r+1}(xi,y)))^(1/(r+1));
            end
        end
        
        v_y = newtonsMethod(@(vy)g(xi,vy)-gdy0-1i*q,@(vy)dgdy{1}(xi,vy),v_y,maxItrs,Eps);
        dv_y = 1i/dgdy{1}(xi,v_y);
        
        F = F + expFact*f(xi,v_y)*dv_y*fact;
    end
end
F = F/omega;

function F = Fyx_(x,y,r,f,g,dgdx,dgdy,dirx,jsy,omega,noGLpts,plotPaths,useFreudQuad,cp,Eps,maxItrs,exactPath,u_x_,v_y_)

if useFreudQuad
    [Q1D_p, W1D_p] = gaussFreudQuad(noGLpts,r(1));
    [Q1D_q, W1D_q] = gaussFreudQuad(noGLpts,r(2));
else
    [Q1D_p, W1D_p] = gaussLaguerreQuad(noGLpts,r(1));
    [Q1D_q, W1D_q] = gaussLaguerreQuad(noGLpts,r(2));
    alpha = -r./(1+r);
end

gdx0 = g(x,y);
if r(1) == 0
    gdx1 = dgdx{1}(x,y);
    gdx2 = dgdx{2}(x,y);
end
F = 0;
for i = 1:size(W1D_p,1)
    wt_p = W1D_p(i);
    if useFreudQuad
        ptt = Q1D_p(i);
        pt = ptt^(r(1)+1);
    else
        pt = Q1D_p(i);
    end
    p = pt/omega;
    if exactPath
        u_x = u_x_(p);
    else
        if r(1) == 0
            if gdx2 == 0
                u_x = x + 1i*p/gdx1;
            else
                u_x = x - (gdx1-sign(real(gdx1))*sqrt(gdx1^2+2*1i*gdx2*p))/gdx2;
            end
        else
            u_x = x + dirx*(factorial(r(1)+1)*p/abs(dgdx{r(1)+1}(x,y)))^(1/(r(1)+1));
        end
        [u_x,noitr] = newtonsMethod(@(ux)g(ux,y)-gdx0-1i*p,@(ux)dgdx{1}(ux,y),u_x,maxItrs,Eps);
    %     noitr
    end
    du_x = 1i/dgdx{1}(u_x,y);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotPaths
        figure(42)
        hold on
        plot(u_x,'*','color','red')
%         keyboard
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gdy0 = g(u_x,y);
    if r(2) == 0
        gdy1 = dgdy{1}(u_x,y);
        gdy2 = dgdy{2}(u_x,y);
%             gdy3 = dgdy{3}(u_x,y);
%             gdy4 = dgdy{4}(u_x,y);
    end
    
    for j = 1:size(W1D_q,1)
        wt = W1D_q(j)*wt_p;
        if useFreudQuad
            qtt = Q1D_q(j);
            qt = qtt^(r(2)+1);
        else
            qt = Q1D_q(j);
        end
        
        q = qt/omega;
        if useFreudQuad
            fact = (r(1)+1)*(r(2)+1)*ptt^r(1)*qtt^r(2)*wt;
        else
            fact = pt^(-alpha(1))*qt^(-alpha(2))*wt;
        end
        
        if exactPath
            v_y = v_y_(u_x,q);
        else
            if r(2) == 0
                if gdy2 == 0
                    v_y = y + 1i*q/gdy1;
                else
                    v_y = y - (gdy1-sign(real(gdy1))*sqrt(gdy1^2+2*1i*gdy2*q))/gdy2;
                end
    %             v_y = y + 1i*q/gdy1 + gdy2/2/gdy1^3*q^2;
    %             v_y = y + 1i*q/gdy1 + gdy2/2/gdy1^3*q^2 + (gdy3*gdy1-3*gdy2^2)/6/gdy1^5*1i*q^3;
    %             v_y = y + 1i*q/gdy1 + gdy2/2/gdy1^3*q^2 + (gdy3*gdy1-3*gdy2^2)/6/gdy1^5*1i*q^3 ...
    %                     - (gdy4*gdy1^2-10*gdy3*gdy2*gdy1+15*gdy2^3)/24/gdy1^7*q^4;

    %             P = 6*gdy1/gdy3-3*(gdy2/gdy3)^2;
    %             Q = 6*gdy1*gdy2/gdy3^2 + 6*1i*q/gdy3-2*(gdy2/gdy3)^3;
    %             jj = 2;
    %             W = exp(2*pi*1i*jj/3)*(Q/2+sqrt(Q^2/4+P^3/27))^(1/3);
    %             v_y2 = y + W - P/(3*W) - gdy2/gdy3; 
            else
                diry = findDirsNSD2(jsy,cp,dgdx,dgdy,2,u_x);
                v_y = y + diry(cp(1,:) == y)*(factorial(r(2)+1)*q/abs(dgdy{r(2)+1}(u_x,y)))^(1/(r(2)+1));
            end

            [v_y,noitr] = newtonsMethod(@(vy)g(u_x,vy)-gdy0-1i*q,@(vy)dgdy{1}(u_x,vy),v_y,maxItrs,Eps);
    %         noitr
        end
        dv_y = 1i/dgdy{1}(u_x,v_y);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if plotPaths
            figure(52)
            plot(v_y,'*','color','red')
            hold on
%             keyboard
%             plotFractal(x,y,u_x,v_y,r,f,g,dgdx,dgdy,q)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        F = F + f(u_x,v_y)*dv_y*du_x*fact;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if plotPaths
%         figure(42)
%         hold on
%         nppts = 1000;
% 
%         x = linspace(0,2,nppts);
%         y = linspace(-6,6,nppts);
%         [X,Y] = meshgrid(x,y);
%         Z = imag(1i*g(u_x,X+1i*Y));
%         lvlList = 10.^linspace(log10(min(min(abs(Z)))),log10(max(max(abs(Z)))),20);
%         lvlList = [-fliplr(lvlList), lvlList];
%         
%         [~,h1] = contour(x,y,Z,40,'LevelList',lvlList);
%         uistack(h1, 'bottom');
%         colorbar
%         drawnow
%         hold off
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotPaths
        nppts = 1000;
        if useFreudQuad
            q_end = Q1D_q(end)^(r(2)+1)/omega;
        else
            q_end = Q1D_q(end)/omega;
        end
        q_plot = linspace(Eps,q_end,nppts);
        if r(2) == 0
            gdy1 = dgdy{1}(u_x,y);
            gdy2 = dgdy{2}(u_x,y);
            if gdy2 == 0
                v_y_plot2 = y + 1i*q_plot/gdy1;
                v_y = y + 1i*q_plot(1)/gdy1;
            else
                v_y_plot2 = y - (gdy1-sign(real(gdy1))*sqrt(gdy1^2+2*1i*gdy2*q_plot))/gdy2;
                v_y = y - (gdy1-sign(real(gdy1))*sqrt(gdy1^2+2*1i*gdy2*q_plot(1)))/gdy2;
            end
        else
            v_y_plot2 = y - dirs*(factorial(r(2)+1)*q_plot/abs(dgdy{r(2)+1}(u_x,y)))^(1/(r(2)+1)); 
            v_y = y + dirs*(factorial(r(2)+1)*q_plot(1)/abs(dgdy{r(2)+1}(u_x,y)))^(1/(r(2)+1));
        end
        v_y_plot = zeros(1,nppts);
        for j = 1:length(q_plot)
            v_y = newtonsMethod(@(vy)g(u_x,vy)-g(u_x,y)-1i*q_plot(j),@(vy)dgdy{1}(u_x,vy),v_y,maxItrs,Eps);
            v_y_plot(j) = v_y;      
        end
        plot(v_y_plot,'red')
        plot(v_y_plot2,'--','color','blue')
%         keyboard
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
F = F*exp(1i*omega*g(x,y))/omega^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotPaths
    hold off
    figure(42)
    hold on
    nppts = 1000;
    if useFreudQuad
        p_end = Q1D_p(end)^(r(1)+1)/omega;
    else
        p_end = Q1D_p(end)/omega;
    end
    p_plot = linspace(Eps,p_end,nppts);
    p = p_plot(end);
    if r(1) == 0
        gdx1 = dgdx{1}(x,y);
        gdx2 = dgdx{2}(x,y);
        if gdx2 == 0
            u_x_plot2 = x + 1i*p/gdx1;
            u_x = x + 1i*p/gdx1;
        else
            u_x_plot2 = x - (gdx1-sign(real(gdx1))*sqrt(gdx1^2+2*1i*gdx2*p_plot))/gdx2;
            u_x = x - (gdx1-sign(real(gdx1))*sqrt(gdx1^2+2*1i*gdx2*p))/gdx2;
        end
    else
        u_x_plot2 = x + dirx*(factorial(r(1)+1)*p_plot/abs(dgdx{r(1)+1}(x,y))).^(1/(r(1)+1));
        u_x       = x + dirx*(factorial(r(1)+1)*p/abs(dgdx{r(1)+1}(x,y)))^(1/(r(1)+1));
    end
    u_x_plot = zeros(1,nppts);
    for j = 1:length(p_plot)
        p = p_plot(end-j+1);
        u_x = newtonsMethod(@(ux)g(ux,y)-g(x,y)-1i*p,@(ux)dgdx{1}(ux,y),u_x,maxItrs,Eps);
        if abs(g(u_x,y)-g(x,y)-1i*p) > Eps
            break
        end
        u_x_plot(j) = u_x;      
    end
    u_x_plot(j:end) = [];
    plot(u_x_plot,'red')
    plot(u_x_plot2,'--','color','blue')
%     keyboard
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if task
%     gstr = func2str(g);
%     nppts = 1000;
%     xLims = xlim;
%     xLims = xLims + [-1,1]*0.1*max(abs(xLims));
%     yLims = ylim;
%     yLims = yLims + [-1,1]*0.1*max(abs(yLims));
%     yLims = [-0.25,0];
%     xLims = [-1.2, 1.2];
%     
%     x = linspace(xLims(1),xLims(2),nppts);
%     y = linspace(yLims(1),yLims(2),nppts);
%     [X,Y] = meshgrid(x,y);
%     Z = imag(1i*g(u_x,X+1i*Y));
% 
% %     if strcmp(gstr(5:end),'x')
%     maxZ = max(max(abs(Z)));   
%     lvlstart = -1.5;
%     lvlList = 10.^linspace(lvlstart,log10(maxZ),40);
%     lvlList = [-fliplr(lvlList), lvlList];
% %         contour(x,y,Z,'LevelList',lvlList)
% 
%     [~,h1] = contour(x,y,Z,40);
%     uistack(h1, 'bottom');
% %     end
%     fstr = func2str(f);
%     title(['g(x) = ' gstr(7:end) ' f(x) = ' fstr(7:end) ', \omega = ' num2str(omega)])
%     if task == 14
%         xlim(xLims)
%         ylim(yLims)
%     end
%     colorbar
%     drawnow
%     hold off
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function plotFractal(x,y,u_x,v_y,r,f,g,dgdx,dgdy,q)
%                      
% figure(62)
% xLims = [-1.2,-0.8];
% yLims = [-0.05,0];
% 
% npxs = 1000;
% % npxs = 20000;
% x_arr = linspace(xLims(1),xLims(2),npxs);
% y_arr = linspace(yLims(1),yLims(2),npxs*(yLims(2)-yLims(1))/(xLims(2)-xLims(1)));
% Eps = 1e-3*(xLims(2)-xLims(1))/npxs;
% 
% max_iteration = 400;
% itrs = zeros(length(y_arr),length(x_arr));
% for ii = 1:length(x_arr)
%     x = x_arr(ii);
%     fprintf('Completed %d out of %d columns!\n', ii, length(x_arr))
%     parfor jj = 1:length(y_arr)
%         z = x + 1i*y_arr(jj);
%         [v_y, iteration] = newtonsMethod(@(vy)g(u_x,vy)-g(u_x,y)-1i*q,@(vy)dgdy(u_x,vy),z,max_iteration,Eps);
%         
%         phi = abs(g(u_x,v_y)-g(u_x,y)-1i*q)/abs(1i*q);
% %         if phi < 10*Eps
%             itrs(jj,ii) = iteration;
% %         else
% %             itrs(jj,ii) = Inf;
% %         end
%     end
% end
% itrs(itrs == Inf) = 401;
% if 1
%     itrs_plot = log10(itrs);
% else
%     itrs_plot = itrs;
% end
% figure(52)
% imagesc(x_arr,y_arr,itrs_plot)
% set(gca,'YDir','normal')
% % setFractalColormap(1000,2,1.7)
% colormap default
% colorbar
% axis equal tight
% hold on
% plot(v_y,'*')
% % savefig('../graphics/largeGraphicsFiles/fractal_g_x3.fig')
