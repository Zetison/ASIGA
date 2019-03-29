function color = uTest_BC(v,a,b,c,L,g2,g3,alpha,beta,s,part,nurbs)


switch part
    case 1 % backpart
        h = g2*tan(alpha/2); % = g2/tan((pi-alpha)/2)

        x2 = g3*tan(alpha);
        r = sqrt(v(:,2).^2+v(:,3).^2);
        r0 = abs(-b+h+x2);
        x0 = -(L+g2+g3);
        indices = r > r0;
        color = zeros(size(v,1),1);
        color(indices) = sqrt((r(indices)-r0).^2+(v(indices,1)-x0).^2);
        color(~indices) = abs(v(~indices,1)-x0);
        return
%     case 2 % cone
%         h = g2*tan(alpha/2); % = g2/tan((pi-alpha)/2)
%         L_c = -(L+g2+cot(alpha)*(b-h));
%         color = abs(norm2(v(:,2:3))*cos(alpha) - (v(:,1)-L_c)*sin(alpha));
%         return
    case 3 % lower transition part
        c1 = g2*csc(alpha);
        c2 = c1-b;
                
        eta = atan2(v(:,3),v(:,2));
        xi = atan2(-cos(eta).*v(:,2)-sin(eta).*v(:,3)-c2, L+v(:,1));
        
        X = @(xi,eta) [-L+c1*cos(xi), -(c2+c1*sin(xi)).*cos(eta), -(c2+c1*sin(xi)).*sin(eta)];
%         % do one steps of Newtons method to reduce round-off errors
%         J{1} = @(xi,eta) -norm2(dXdxi(xi,eta)).^2 + sum((v - X(xi,eta)).*d2Xdxi2(xi,eta),2);
%         J{2} = @(xi,eta) -sum(dXdxi(xi,eta).*dXdeta(xi,eta),2) + sum((v - X(xi,eta)).*d2Xdxideta(xi,eta),2);
%         J{3} = @(xi,eta) -norm2(dXdeta(xi,eta)).^2 + sum((v - X(xi,eta)).*d2Xdeta2(xi,eta),2);
% 
%         F{1} = @(xi,eta) sum((v - X(xi,eta)).*dXdxi(xi,eta),2);
%         F{2} = @(xi,eta) sum((v - X(xi,eta)).*dXdeta(xi,eta),2);
%         [xi,eta] = newtonsMethod2D(F,J,[xi,eta],1,1e-15);
        color = norm2(v-X(xi,eta));
    case num2cell(4:15) % upper transition part
%         color = zeros(size(v,1),1);
%         for i = 1:numel(color)
%             color(i) = minDistNURBS(v(i,:).',nurbs);
%         end
        FJ = @(xi,eta) evalNURBSnewton(xi,eta,v,nurbs);
        if v(round(end/2),1) > -45.290511658811013
            eta = 3*ones(size(v,1),1)/4;
            bnd = [0,1;
                   0.5,1];
        else
            eta = ones(size(v,1),1)/4;
            bnd = [0,1;
                   0,0.5];
        end
        xi = ones(size(v,1),1)/2;
        [xi,eta] = newtonsMethod2D(FJ,NaN,[xi,eta],10,1e-15,bnd);
        X = evaluateNURBS_2ndDeriv(nurbs, [xi,eta]);
        color = norm2(v-X);
    case 16 %{127,189} % lower cylindrical part
        color = abs(norm2(v(:,2:3)) - b);
    case 22 %num2cell([150:153,163:166]) % deck
        color = abs(v(:,3) - c);
    case num2cell([17:21,23:27]) % side panels
        C_4 = c + b*cos(beta/2);
        C_3 = b*sin(beta/2)-s;
        C_2 = (2*C_4+C_3*tan(beta/2))/C_3^3;
        C_1 = -(3*C_4+C_3*tan(beta/2))/C_3^2;

        npts = 6;
        dy = C_3/(npts-1);
        P = @(y) c + C_1*(y-s).^2 + C_2*(y-s).^3;
        if part > 22
            sgn = -1;
        else
            sgn = 1;
        end
        if part == 17 || part == 27 %any(ismember(part,[128,188])) % 21, 31
            i = 4;
        elseif part == 18 || part == 26 %any(ismember(part,[129,187])) % 22, 30
            i = 3;
        elseif part == 19 || part == 25 %any(ismember(part,[130,186])) % 23, 29
            i = 2;
        elseif part == 20 || part == 24 %any(ismember(part,[131,133:135,148,168,181:183,185])) % 24, 28
            i = 1;
        else  % 25, 27
            i = 0;
        end
        y0 = s+i*dy;
        y1 = s+(i+1)*dy;
        z0 = P(y0);
        z1 = P(y1);
        P1 = [0, sgn*y0, z0];
        P2 = [0, sgn*y1, z1];
        P3 = [-L, sgn*y1, z1];
        normal = sgn*cross(P3-P2,P1-P2);
        D1 = normal(2);
        D2 = normal(3);
        D3 = D1*P1(2) + D2*P1(3);
        y = v(:,2);
        z = v(:,3);
        color = abs(D1*y+D2*z-D3)/norm(normal);
    case 28 % bow lower part
        r_yz = sqrt(v(:,2).^2+v(:,3).^2);
        indices = r_yz < 10*eps;
        color = zeros(size(v,1),1);   
        if any(indices)
            color(indices) = abs(a-v(indices,1));
        end
        
        indices = abs(v(:,1)) < 10*eps;   
        if any(indices)
            color(indices) = abs(b-r_yz(indices));
        end
        
        indices = and(r_yz >= 10*eps, abs(v(:,1)) >= 10*eps);        
        if any(indices)
            r_yz = r_yz(indices);
            v = v(indices,:);
            eta = atan2(v(:,3),v(:,2));
            at = -(a^2-b^2)^2;
            bt = -2*(a^2-b^2)*b*r_yz;
            ct = (a^2-b^2)^2-(b^2*r_yz.^2 + a^2*v(:,1).^2);
            dt = 2*(a^2-b^2)*b*r_yz;
            et = b^2*r_yz.^2;
            
            Delta0 = ct.^2 - 3*bt.*dt+12*at*et;
            Delta1 = 2*ct.^3-9*bt.*ct.*dt+27*bt.^2.*et+27*at*dt.^2-72*at*ct.*et;
%             Delta = -(Delta1^2-4*Delta0^3)/27;
            phi = acos(Delta1./(2*sqrt(Delta0.^3)));
            p = (8*at*ct-3*bt.^2)/(8*at^2);
            q = (bt.^3-4*at*bt.*ct+8*at.^2.*dt)/(8*at^3);
            S = sqrt(-2*p/3+2/(3*at).*sqrt(Delta0).*cos(phi/3))/2;
%             sinxi = -bt/(4*at) - S + [-1,1]*sqrt(-4*S.^2-2*p+q./S)/2;
%             sinxi = -bt/(4*at) + S + [-1,1]*sqrt(-4*S.^2-2*p-q./S)/2;
            
            sinxi = -bt/(4*at) + S + sqrt(-4*S.^2-2*p-q./S)/2;
            xi = abs(asin(sinxi));
            X = @(xi,eta) [a*cos(xi), b*sin(xi).*cos(eta), b*sin(xi).*sin(eta)];
            dXdxi = @(xi,eta) [-a*sin(xi), b*cos(xi).*cos(eta), b*cos(xi).*sin(eta)];
            dXdeta = @(xi,eta) [zeros(size(xi)), -b*sin(xi).*sin(eta), b*sin(xi).*cos(eta)];
            d2Xdxi2 = @(xi,eta) [-a*cos(xi), -b*sin(xi).*cos(eta), b*sin(xi).*sin(eta)];
            d2Xdeta2 = @(xi,eta) [zeros(size(xi)), -b*sin(xi).*cos(eta), -b*sin(xi).*sin(eta)];
            d2Xdxideta = @(xi,eta) [zeros(size(xi)), -b*cos(xi).*sin(eta), b*cos(xi).*cos(eta)];
            
            % do one steps of Newtons method to reduce round-off errors
            J{1} = @(xi,eta) -norm2(dXdxi(xi,eta)).^2 + sum((v - X(xi,eta)).*d2Xdxi2(xi,eta),2);
            J{2} = @(xi,eta) -sum(dXdxi(xi,eta).*dXdeta(xi,eta),2) + sum((v - X(xi,eta)).*d2Xdxideta(xi,eta),2);
            J{3} = @(xi,eta) -norm2(dXdeta(xi,eta)).^2 + sum((v - X(xi,eta)).*d2Xdeta2(xi,eta),2);

            F{1} = @(xi,eta) sum((v - X(xi,eta)).*dXdxi(xi,eta),2);
            F{2} = @(xi,eta) sum((v - X(xi,eta)).*dXdeta(xi,eta),2);
            [xi,eta] = newtonsMethod2D(F,J,[xi,eta],1,1e-15);

            color(indices) = norm2(v-X(xi,eta));
        end
    case num2cell(29:39)  % upper bow
        C_4 = c + b*cos(beta/2);
        C_3 = b*sin(beta/2)-s;
        C_2 = (2*C_4+C_3*tan(beta/2))/C_3^3;
        C_1 = -(3*C_4+C_3*tan(beta/2))/C_3^2;

        npts = 6;
        dy = C_3/(npts-1);
        P = @(y) c + C_1*(y-s).^2 + C_2*(y-s).^3;
        
        
        if part > 34
            sgn = -1;
            i = part-35;
        else
            sgn = 1;
            i = 33-part;
        end
        if part == 34
            y1 = s;
            y2 = -s;
            z1 = c;
            z2 = c;
        else
            y1 = sgn*(s+i*dy);
            y2 = sgn*(s+(i+1)*dy);
            z1 = P(abs(y1));
            z2 = P(abs(y2));
        end
        theta1 = atan2(z1,y1);
        theta2 = atan2(z2,y2);
        r1 = norm([y1,z1]);
        r2 = norm([y2,z2]);
        
        C1 = r2*cos(theta2);
        C2 = r1*cos(theta1)-r2*cos(theta2);
        C3 = r2*sin(theta2);
        C4 = r1*sin(theta1)-r2*sin(theta2);
        
        r_yz = sqrt(v(:,2).^2+v(:,3).^2);
        indices = r_yz < 10*eps;
        color = zeros(size(v,1),1);   
        if any(indices)
            color(indices) = abs(a-v(indices,1));
        end
        X = @(xi,eta) [a*cos(xi), (C1 + eta*C2).*sin(xi), (C3 + eta*C4).*sin(xi)];

        indices = abs(v(:,1)) < 10*eps;   
        if any(indices)
            eta = -(C1*C2-C2*v(indices,2)+C3*C4-C4*v(indices,3))/(C2^2+C4^2);
            color(indices) = norm2(v(indices,:)-X(pi/2*ones(numel(eta),1),eta));
        end

        indices = and(abs(v(:,2)) < 10*eps,~or(r_yz < 10*eps, abs(v(:,1)) < 10*eps));   
        if any(indices)
            v1 = v(indices,1);
            v3 = v(indices,3);
            xi = asin(abcdeFormula([(C3^2-a^2)^2*ones(size(v1)), 2*C3*v3*(a^2-C3^2), C3^2*v3.^2-(a^2-C3^2)^2+a^2*v1.^2, ...
                               -2*C3*v3*(a^2-C3^2), -C3^2*v3.^2]));
            eta = 0.5*ones(size(xi));

            temp = zeros(size(v1,1),4);
            for i = 1:4
                temp(:,i) = norm2(v(indices,:)-X(xi(:,i),eta(:,i)));
            end
            xiTemp = xi;
            etaTemp = eta;
            xi = zeros(size(v1,1),1);
            eta = zeros(size(v1,1),1);
            for j = 1:size(v1,1)
                [~,I] = min(temp(j,:));
                xi(j) = xiTemp(j,I);
                eta(j) = etaTemp(j,I);
            end

            color(indices) = norm2(v(indices,:)-X(xi,eta));
        end

        indices = and(and(r_yz >= 10*eps, abs(v(:,1)) >= 10*eps),abs(v(:,2)) >= 10*eps); 
%         test = C4*C3+C1*C2;
        v_old = v;
        if any(indices)
            v = v(indices,:);

            C5 = -(C2*v(:,3)-C4*v(:,2))*(C1*C4-C2*C3);
            C6 = -C1^2*C4*v(:,3)+C3*(C2*v(:,3)+C4*v(:,2))*C1+(-C3^2*v(:,2)+a^2*v(:,2))*C2+C4*a^2*v(:,3);
            C7 = (C2^2+C4^2)^2;
            C8 = (2*(C1*C2+C3*C4))*(C2^2+C4^2);
            C9 = ((C1-v(:,2))*C2+C4*(C3-v(:,3))).*((C1+v(:,2))*C2+C4*(C3+v(:,3)));
            C10 = (C2*v(:,2)+C4*v(:,3)).^2*a^2.*v(:,1).^2*(C2^2+C4^2)^2;
            C11 = 2*(C2*v(:,2)+C4*v(:,3)).^2*a^2.*v(:,1).^2*(C1*C2+C3*C4)*(C2^2+C4^2);
            C12 = (C2*v(:,2)+C4*v(:,3)).^2*a^2.*v(:,1).^2*(C1*C2+C3*C4)^2;

            eta = abcdeFormula([C5.^2*C7, C5.^2*C8+2*C5.*C6*C7, C5.^2.*C9+2*C5.*C6*C8+C6.^2*C7-C10, 2*C5.*C6.*C9+C6.^2*C8-C11, C6.^2.*C9-C12]);
            xi = asin((C2*v(:,[2,2,2,2])+C4*v(:,[3,3,3,3]))./(C2^2*eta+C1*C2+C4*(C4*eta+C3)));
            indices2 = and(and(-10*eps<=xi,xi<=pi/2+10*eps),and(-10*eps<=eta,eta<=1+10*eps));
%                 xi = xi(:,4);
%                 eta = eta(:,4);
%                 eta = zeros(size(etaTemp,1),1);
%                 indices2 = and(0 <= etaTemp(:,1), etaTemp(:,1) <= 1);
%                 eta(indices2) = etaTemp(indices2,1);
%                 eta(~indices2) = etaTemp(~indices2,2);

            X = @(xi,eta) [a*cos(xi), (C1 + eta*C2).*sin(xi), (C3 + eta*C4).*sin(xi)];
            dXdxi = @(xi,eta) [-a*sin(xi), (C1 + eta*C2).*cos(xi), (C3 + eta*C4).*cos(xi)];
            dXdeta = @(xi,eta) [zeros(size(xi)), C2.*sin(xi), C4.*sin(xi)];
            d2Xdxi2 = @(xi,eta) [-a*cos(xi), -(C1 + eta*C2).*sin(xi), -(C3 + eta*C4).*sin(xi)];
            d2Xdeta2 = @(xi,eta) [zeros(size(xi)), zeros(size(xi)), zeros(size(xi))];
            d2Xdxideta = @(xi,eta) [zeros(size(xi)), C2.*cos(xi), C4.*cos(xi)];

            temp = zeros(size(v,1),4);
            for i = 1:4
                temp(:,i) = norm2(v-X(xi(:,i),eta(:,i)));
            end
            xiTemp = xi;
            etaTemp = eta;
            xi = zeros(size(v,1),1);
            eta = zeros(size(v,1),1);
            for j = 1:size(v,1)
                [~,I] = min(temp(j,:));
                xi(j) = xiTemp(j,I);
                eta(j) = etaTemp(j,I);
            end
            % do one step of Newtons method to reduce round-off errors
            J{1} = @(xi,eta) -norm2(dXdxi(xi,eta)).^2 + sum((v - X(xi,eta)).*d2Xdxi2(xi,eta),2);
            J{2} = @(xi,eta) -sum(dXdxi(xi,eta).*dXdeta(xi,eta),2) + sum((v - X(xi,eta)).*d2Xdxideta(xi,eta),2);
            J{3} = @(xi,eta) -norm2(dXdeta(xi,eta)).^2 + sum((v - X(xi,eta)).*d2Xdeta2(xi,eta),2);

            F{1} = @(xi,eta) sum((v - X(xi,eta)).*dXdxi(xi,eta),2);
            F{2} = @(xi,eta) sum((v - X(xi,eta)).*dXdeta(xi,eta),2);
            [xi,eta] = newtonsMethod2D(F,J,[xi,eta],5,1e-15);
            color(indices) = norm2(v-X(xi,eta));
        end
%             XI = linspace(0,pi/2,100);
%             ETA = linspace(0,1,100);
%             [XI,ETA] = meshgrid(XI,ETA);
%             X = a*cos(XI);
%             Y = (C1 + ETA*C2).*sin(XI);
%             Z = (C3 + ETA*C4).*sin(XI);
% %             surf(X,Y,Z, 'EdgeColor','none','LineStyle','none')
%             surf(X,Y,Z)

        X = @(t) [a*cos(t(1,:)); (C1 + t(2,:)*C2).*sin(t(1,:));
                                 (C3 + t(2,:)*C4).*sin(t(1,:))];
        dXdxi = @(t) [-a*sin(t(1,:)); (C1 + t(2,:)*C2).*cos(t(1,:));
                                      (C3 + t(2,:)*C4).*cos(t(1,:))];
        dXdeta = @(t) [zeros(1,size(t,2)); C2*sin(t(1,:));
                                           C4*sin(t(1,:))];
        d2Xdxi2 = @(t) [-a*cos(t(1,:)); -(C1 + t(2,:)*C2).*sin(t(1,:));
                                        -(C3 + t(2,:)*C4).*sin(t(1,:))];
        d2Xdeta2 = @(t) zeros(3,size(t,2));
        d2Xdxideta = @(t) [zeros(1,size(t,2)); C2*cos(t(1,:));
                                               C4*cos(t(1,:))];
        bnd = [0,pi/2;
               0,1];
        indices = color > 1e-14;
        v = v_old(indices,:);
        colorTemp = zeros(size(v,1),1);
        parfor i = 1:size(v,1)
            colorTemp(i) = minDist(v(i,:).',X,dXdxi,dXdeta,d2Xdxi2,d2Xdeta2,d2Xdxideta,bnd);
        end
        color(indices) = colorTemp;
    case 42 %{157,158,159} % sail roof
        h_s = 3.5;
        color = abs(v(:,3) - (c+h_s));
    case {41,43} % sail
        R_x = rotationMatrix(0, 'Xaxis');
        if part == 43
            sgn = -1;
        else
            sgn = 1;
        end
        b_ls = 2*s;
        b_us = 2;
        l_ls = 13;
        l_us = 12.3;
        h_s = 3.5;
        delta_s = 0.2;
        T = [a-19,0,c].';
        [f_u,dfxi_u,dfdxi2_u] = NACA(b_ls,l_ls);
        [f_o,dfxi_o,dfdxi2_o] = NACA(b_us,l_us);
        [X,dXdxi,dXdeta,d2Xdxi2,d2Xdeta2,d2Xdxideta] = rudder(l_ls,l_us,h_s,delta_s,R_x,T,sgn,f_u,dfxi_u,dfdxi2_u,f_o,dfxi_o,dfdxi2_o);
        bnd = [0,1;
               0,1];
        FJ = @(xi,eta) evalNewton(xi,eta,v,X,dXdxi,dXdeta,d2Xdxi2,d2Xdeta2,d2Xdxideta);
        xi = ones(size(v,1),1)*0.9;
        eta = ones(size(v,1),1)*mean(bnd(2,:));
        [xi,eta] = newtonsMethod2D(FJ,NaN,[xi,eta],15,1e-15,bnd);
        color = norm2(v-X([xi,eta].').');
%         color = zeros(size(v,1),1);
%         parfor i = 1:numel(color)
%             color(i) = minDist(v(i,:).',X,dXdxi,dXdeta,d2Xdxi2,d2Xdeta2,d2Xdxideta,bnd);
%         end
    case 45 %num2cell(17:21) % rudder top
        h_s = 3.5;
        color = abs(v(:,3) - h_s);
    case 48 %num2cell(43:47) % rudder righ_s
        h_s = 3.5;
        color = abs(-v(:,2) - h_s);
    case 51 %num2cell(69:73) % rudder bottom
        h_s = 3.5;
        color = abs(-v(:,3) - h_s);
    case 54 %num2cell(95:99) % rudder left
        h_s = 3.5;
        color = abs(v(:,2) - h_s);
    case {2,44,46,47,49,50,52,53,55} % main rudders
        if part == 44 || part == 46
            R_x = rotationMatrix(0, 'Xaxis');
            if part == 46
                sgn = -1;
            else
                sgn = 1;
            end
        elseif part == 47 || part == 49
            R_x = rotationMatrix(pi/2, 'Xaxis');
            if part == 49
                sgn = -1;
            else
                sgn = 1;
            end
        elseif part == 50 || part == 52
            R_x = rotationMatrix(pi, 'Xaxis');
            if part == 52
                sgn = -1;
            else
                sgn = 1;
            end
        else
            R_x = rotationMatrix(-pi/2, 'Xaxis');
            if part == 55
                sgn = -1;
            else
                sgn = 1;
            end
        end
        l_ls = 2.6;
        b_ls = 0.4;
        l_us = 2.35;
        b_us = 0.3;
        h_s = 3.5;
        delta_s = 0.25;
        T = [a-58.9,0,0].';
        [f_u,dfxi_u,dfdxi2_u] = NACA(b_ls,l_ls);
        [f_o,dfxi_o,dfdxi2_o] = NACA(b_us,l_us);
        [X,dXdxi,dXdeta,d2Xdxi2,d2Xdeta2,d2Xdxideta] = rudder(l_ls,l_us,h_s,delta_s,R_x,T,sgn,f_u,dfxi_u,dfdxi2_u,f_o,dfxi_o,dfdxi2_o);
        bnd = [0,1;
               0,1];
        FJ = @(xi,eta) evalNewton(xi,eta,v,X,dXdxi,dXdeta,d2Xdxi2,d2Xdeta2,d2Xdxideta);
        xi = ones(size(v,1),1)*0.9;
        eta = ones(size(v,1),1)*mean(bnd(2,:));
        [xi,eta] = newtonsMethod2D(FJ,NaN,[xi,eta],15,1e-15,bnd);
        color = norm2(v-X([xi,eta].').');
        
        
        h = g2*tan(alpha/2); % = g2/tan((pi-alpha)/2)
        L_c = -(L+g2+cot(alpha)*(b-h));
        color = min([color, abs(norm2(v(:,2:3))*cos(alpha) - (v(:,1)-L_c)*sin(alpha))],[],2);
        return
%         color = zeros(size(v,1),1);
%         parfor i = 1:numel(color)
%             color(i) = minDist(v(i,:).',X,dXdxi,dXdeta,d2Xdxi2,d2Xdeta2,d2Xdxideta,bnd);
%         end
    case 57 %{139,140,141} % depth rudder left
        color = abs(v(:,2) - b);
    case 60 % depth rudder right
        color = abs(-v(:,2) - b);
    case {56,58,59,61} % depth rudders
        if part == 56 || part == 58 %ismember(part,[136:138,142:144])
            R_x = rotationMatrix(-pi/2, 'Xaxis');
            if part == 58 %140
                sgn = -1;
            else
                sgn = 1;
            end
            sgn2 = 1;
        else
            R_x = rotationMatrix(pi/2, 'Xaxis');
            if part == 61 %176
                sgn = -1;
            else
                sgn = 1;
            end
            sgn2 = -1;
        end
        C_4 = c + b*cos(beta/2);
        C_3 = b*sin(beta/2)-s;
        C_2 = (2*C_4+C_3*tan(beta/2))/C_3^3;
        C_1 = -(3*C_4+C_3*tan(beta/2))/C_3^2;

        npts = 6;
        dy = C_3/(npts-1);
        P = @(y) c + C_1*(y-s).^2 + C_2*(y-s).^3;

        l_ls = 2.6;
        b_ls = 2*(c-P(s+dy));
        l_us = 2.35;
        b_us = 0.22;
        h_s = b-s;
        delta_s = 0.25;
        T = [a-11.0,sgn2*s,c-b_ls/2].';
        [f_u,dfxi_u,dfdxi2_u] = NACA(b_ls,l_ls);
        [f_o,dfxi_o,dfdxi2_o] = NACA(b_us,l_us);
        [X,dXdxi,dXdeta,d2Xdxi2,d2Xdeta2,d2Xdxideta] = rudder(l_ls,l_us,h_s,delta_s,R_x,T,sgn,f_u,dfxi_u,dfdxi2_u,f_o,dfxi_o,dfdxi2_o);
        bnd = [0,1;
               0,1];
        FJ = @(xi,eta) evalNewton(xi,eta,v,X,dXdxi,dXdeta,d2Xdxi2,d2Xdeta2,d2Xdxideta);
        xi = ones(size(v,1),1)*0.9;
        eta = ones(size(v,1),1)*mean(bnd(2,:));
        [xi,eta] = newtonsMethod2D(FJ,NaN,[xi,eta],15,1e-15,bnd);
        color = norm2(v-X([xi,eta].').');
%         color = zeros(size(v,1),1);
%         parfor i = 1:numel(color)
%             color(i) = minDist(v(i,:).',X,dXdxi,dXdeta,d2Xdxi2,d2Xdeta2,d2Xdxideta,bnd);
%         end
    otherwise
        color = 1e7;
end

function FJ = evalNewton(xi,eta,v,X_,dXdxi_,dXdeta_,d2Xdxi2_,d2Xdeta2_,d2Xdxideta_)
X = X_([xi,eta].').';
dXdxi = dXdxi_([xi,eta].').';
dXdeta = dXdeta_([xi,eta].').';
d2Xdxi2 = d2Xdxi2_([xi,eta].').';
d2Xdeta2 = d2Xdeta2_([xi,eta].').';
d2Xdxideta = d2Xdxideta_([xi,eta].').';

FJ{1} = sum((v - X).*dXdxi,2);
FJ{2} = sum((v - X).*dXdeta,2);

FJ{3} = -norm2(dXdxi).^2 + sum((v - X).*d2Xdxi2,2);
FJ{4} = -sum(dXdxi.*dXdeta,2) + sum((v - X).*d2Xdxideta,2);
FJ{5} = -norm2(dXdeta).^2 + sum((v - X).*d2Xdeta2,2);


function FJ = evalNURBSnewton(xi,eta,v,nurbs)

[X,dXdxi,dXdeta,d2Xdxi2,d2Xdeta2,d2Xdxideta] = evaluateNURBS_2ndDeriv(nurbs, [xi,eta]);

FJ{1} = sum((v - X).*dXdxi,2);
FJ{2} = sum((v - X).*dXdeta,2);

FJ{3} = -norm2(dXdxi).^2 + sum((v - X).*d2Xdxi2,2);
FJ{4} = -sum(dXdxi.*dXdeta,2) + sum((v - X).*d2Xdxideta,2);
FJ{5} = -norm2(dXdeta).^2 + sum((v - X).*d2Xdeta2,2);

