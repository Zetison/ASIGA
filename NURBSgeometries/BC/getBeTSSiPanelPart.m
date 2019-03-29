function SS = getBeTSSiPanelPart(ss,eta,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,xi1,s1,s2,part,pathPart)
if nargin < 18
    pathPart = part;
end
xi = invertNACA2_depth(ss,b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,pathPart,s2);

C_4 = c + b*cos(beta/2);
C_3 = b*sin(beta/2)-s;
C_2 = (2*C_4+C_3*tan(beta/2))/C_3^3;
C_1 = -(3*C_4+C_3*tan(beta/2))/C_3^2;

nPanels = 5;
dy = C_3/nPanels;
P = @(y) c + C_1*(y-s).^2 + C_2*(y-s).^3;

y_f = s+dy;

[g,~,S,~,~,D1,D2,D3] = getDepthRudderFunctions(b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,part(end-4:end));


r1 = @(xi) S(xi,g(xi));

switch part
    case {'upper','lower'}
        indices = xi<=xi1;
        SS = zeros(3,numel(xi));
        if ~(numel(indices) == 1 && ~indices)
            if strcmp(part,'upper')
                Y = @(ss) y_f - ss/s1*dy;
            else
                Y = @(ss) y_f + ss/s1*dy;
            end
            r2 = @(ss) [x_d*ones(1,numel(ss));
                        Y(ss);
                        (D3-D1*Y(ss))/D2];
            X0 = r1(xi(indices));
            X1 = r2(ss(indices));
            SS(:,indices) =  X0 + (X1-X0).*repmat(eta(indices),3,1);
        end
        indices = xi>xi1;
        if ~(numel(indices) == 1 && ~indices)
            X = @(ss) X_(ss,s1,s2,x_d,l_ld);
            
            
            if strcmp(part,'upper')
                r2 = @(ss) [X(ss); s*ones(size(ss)); c*ones(size(ss))];
            else
                r2 = @(ss) [X(ss); (y_f+dy)*ones(size(ss)); P(y_f+dy)*ones(size(ss))];
            end
            
            X0 = r1(xi(indices));
            X1 = r2(ss(indices));
            SS(:,indices) =  X0 + (X1-X0).*repmat(eta(indices),3,1);
        end
    case {'rudderupper','rudderlower'}
        SS = S(xi,eta);
    case {'curveupper','curvelower'}
        SS = r1(xi);
end

function x = X_(ss,s1,s2,x_s,l_ls)
% X = @(s) x_f + (s-s1)/(1-s1)*(x_b-x_f)
x = zeros(size(ss)); 
indices = ss <= s2;
if ~(numel(indices) == 1 && ~indices)
    x_s2 = x_s - s2*l_ls;
    x(indices) = x_s + (ss(indices)-s1)/(s2-s1)*(x_s2-x_s);
end
% indices = and(ss > s2,ss < 0.3);
% if ~(numel(indices) == 1 && ~indices)
%     x(indices) = x_s - ss(indices)*(l_ls-lt_m);
% end
indices = ss > s2;
if ~(numel(indices) == 1 && ~indices)
    x(indices) = x_s - ss(indices)*l_ls;
end