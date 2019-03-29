function SS = getBeTSSiConePart(ss,eta,b,g2,g3,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m,xi1,s1,s2,part,rudder_i)

h = g2*tan(alpha./2); % = g2./tan((pi-alpha)./2)
L_c = L+g2+(b-h).*cot(alpha);

[g,dg,S,dSdxi,dSdeta] = getMainRudderFunctions(b,g2,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m);

xi = invertNACA2(ss,l_lm,delta_m,x_m,g,dg,dSdxi,dSdeta,s2);
r1 = @(xi) S(xi,g(xi));

x_f = -L-g2-g3/2;
theta = 30*pi/180;


switch part
    case 'cone'
        SS = zeros(3,numel(xi));
        indices = xi<=xi1;
        if any(indices)
            if 1
                if rudder_i == 1
                    Xi = [0,0,0,1,1,2,2,3,3,3]/3;
                else
                    Xi = [0,0,0,1,1,1];
                end
                n = numel(Xi)-(2+1);
                controlPtsTemp = parmArc(Xi,30*pi/180);
                controlPts = [x_f*ones(1,n); (x_f+L_c)*tan(alpha)*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
                R_x = rotationMatrix(60*pi/180, 'Xaxis');
                controlPts(1:3,:) = R_x*controlPts(1:3,:);
                nurbs = createNURBSobject(controlPts, Xi);

                r2 = @(s) evaluateNURBS(nurbs,abs((s1-s)/s1));
            else
                Theta = @(s) theta*s/s1;
                r2 = @(s) [x_f*ones(1,numel(s));
                            (x_f+L_c)*tan(alpha)*cos(pi/2-Theta(s));
                            (x_f+L_c)*tan(alpha)*sin(pi/2-Theta(s))];
            end
            X0 = r1(xi(indices));
            X1 = r2(ss(indices));
            x = X0(1,:) + (X1(1,:)-X0(1,:)).*eta(indices);
            y = X0(2,:) + (X1(2,:)-X0(2,:)).*eta(indices);
            z = sqrt(abs((x+L_c).^2*tan(alpha)^2-y.^2));
            SS(:,indices) = [x;y;z];
        end
        indices = and(xi>xi1,ss<s2);
        if any(indices)
            lt_m = g(0)*delta_m;
            X = @(s) X_(s,s1,s2,x_f,l_lm,lt_m,x_m);
            
            
            r2 = @(s) [X(s);
                        (X(s)+L_c)*tan(alpha)*cos(pi/2-theta);
                        (X(s)+L_c)*tan(alpha)*sin(pi/2-theta)];


            X0 = r1(xi(indices));
            X1 = r2(ss(indices));
            
            x = X0(1,:) + (X1(1,:)-X0(1,:)).*eta(indices);
            y = X0(2,:) + (X1(2,:)-X0(2,:)).*eta(indices);
            z = sqrt(abs((x+L_c).^2*tan(alpha)^2-y.^2));
            SS(:,indices) = [x;y;z];
        end
        indices = ss>=s2;
        if any(indices)
            lt_m = g(0)*delta_m;
            X = @(s) X_(s,s1,s2,x_f,l_lm,lt_m,x_m);
            
            
            r2 = @(s) [X(s);
                        (X(s)+L_c)*tan(alpha)*cos(pi/2-theta);
                        (X(s)+L_c)*tan(alpha)*sin(pi/2-theta)];


            X0 = r1(xi(indices));
            X1 = r2(ss(indices));
            npts = size(X0,2);
            Xi = [0,0,0,1,1,1];
            n = numel(Xi)-(2+1);
            etaTemp = eta(indices);
            stemp = ss(indices);
            x = X0(1,:) + (X1(1,:)-X0(1,:)).*etaTemp;
            y = X0(2,:) + (X1(2,:)-X0(2,:)).*eta(indices);
            z = sqrt(abs((x+L_c).^2*tan(alpha)^2-y.^2));
            indices2 = find(indices);
            for i = 1:npts
                x0 = X0(1,i);
%                 if abs(x0 - X1(1,i)) > 100*eps || etaTemp(i) > 0.5
%                     keyboard
%                 end
                theta0 = atan2(X0(3,i),X0(2,i));
                theta1 = atan2(X1(3,i),X1(2,i));
                controlPtsTemp = parmArc(Xi,theta0-theta1);
                controlPts = [x0*ones(1,n); (x0+L_c)*tan(alpha)*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
                R_x = rotationMatrix(theta1, 'Xaxis');
                controlPts(1:3,:) = R_x*controlPts(1:3,:);
                nurbs = createNURBSobject(controlPts, Xi);

                XX = [x(i);y(i);z(i)];
                thetaE1 = atan2(XX(3),XX(2));
                XX = evaluateNURBS(nurbs,1-etaTemp(i));
                thetaE2 = atan2(XX(3),XX(2));
                thetaA = (1-stemp(i))/(1-s2)*thetaE1 + (stemp(i)-s2)/(1-s2)*thetaE2;
%                 isoCurve = evaluateNURBS(nurbs,etaTemp(i));
                SS(:,indices2(i)) = [x0;
                                    (x0+L_c)*tan(alpha)*cos(thetaA);
                                    (x0+L_c)*tan(alpha)*sin(thetaA)];
%                 h = g2*tan(alpha/2); % = g2/tan((pi-alpha)/2)
%                 Sx = -(L+g2+cot(alpha)*(b-h));
%                 XX = SS(:,indices2(i));
%                 if max(abs(norm2(XX(2:3,:).')*cos(alpha) - (XX(1,:).'-Sx)*sin(alpha))) > 1e2*eps
%                     keyboard
%                 end
            end
        end
    case 'rudder'
        SS = S(xi,eta);
    case 'curve'
        SS = r1(xi);
    case 'curve2'
        r2 = @(xi) [x_f*ones(1,numel(xi));
                    (x_f+L_c)*tan(alpha)*cos(pi/2-theta);
                    (x_f+L_c)*tan(alpha)*sin(pi/2-theta)];
        X0 = r1(xi1);
        X1 = r2(xi1);
        x = X0(1,:) + (X1(1,:)-X0(1,:)).*eta;
        y = X0(2,:)+(x-X0(1,:)).*(X1(2,:)-X0(2,:))./(X1(1,:)-X0(1,:));
        z = sqrt(abs((x+L_c).^2*tan(alpha)^2-y.^2));
        SS = [x;y;z];
%     case 'crossCurve'
%         indices = xi<=xi_t;
%         if ~(numel(indices) == 1 && ~indices)
%             if 1
%                 Xi = [0,0,0,1,1,2,2,3,3,3]/3;
%                 controlPtsTemp = parmArc(Xi,30*pi/180);
%                 controlPts = [x_f*ones(1,7); (x_f+L_c)*tan(alpha)*controlPtsTemp(1:2,:); controlPtsTemp(3,:)];
%                 R_x = rotationMatrix(60*pi/180, 'Xaxis');
%                 controlPts(1:3,:) = R_x*controlPts(1:3,:);
%                 nurbs = createNURBSobject(controlPts, Xi);
% 
%                 r2 = @(xi) evaluateNURBS(nurbs,(xi_t-xi)/xi_t);
%             else
%                 Theta = @(xi) theta*xi/xi_t;
%                 r2 = @(xi) [x_f*ones(1,numel(xi));
%                             (x_f+L_c)*tan(alpha)*cos(pi/2-Theta(xi));
%                             (x_f+L_c)*tan(alpha)*sin(pi/2-Theta(xi))];
%             end
%             X0 = r1(xi(indices));
%             X1 = r2(xi(indices));
%             x = X0(1,:) + (X1(1,:)-X0(1,:)).*xi;
%             y = X0(2,:)+(x-X0(1,:)).*(X1(2,:)-X0(2,:))./(X1(1,:)-X0(1,:));
%             z = sqrt(abs((x+L_c).^2*tan(alpha)^2-y.^2));
%             SS = zeros(3,numel(xi));
%             SS(:,indices) = [x;y;z];
%         end
end

function x = X_(s,s1,s2,x_f,l_ls,lt_m,x_s)
% X = @(s) x_f + (s-s1)/(1-s1)*(x_b-x_f)
x = zeros(size(s));
indices = s > s2;
if ~(numel(indices) == 1 && ~indices)
    x(indices) = x_s-lt_m - s(indices)*(l_ls-lt_m);
end
indices = s <= s2;
if ~(numel(indices) == 1 && ~indices)
    x_s2 = x_s-lt_m - s2*(l_ls-lt_m);
    x(indices) = x_f + (s(indices)-s1)/(s2-s1)*(x_s2-x_f);
end