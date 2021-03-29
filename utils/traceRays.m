function task = traceRays(task)

d_inc = task.d_vec.';
O = zeros(size(task.varCol{1}.o,1),3,2);

objectHit = zeros(size(task.varCol{1}.o,1),1,'logical');

N = size(task.varCol{1}.o,1);
O(:,:,1) = task.varCol{1}.o;

% plotRays = 1;
switch task.misc.model
    case 'S1'
        R = task.varCol{1}.R_i;
        d = zeros(N,3);
        for n = 1:size(O,1)
        % parfor n = 1:size(O,1)
            Otemp = O(n,:,:);
            o = Otemp(1,:,1);
            do = dot(d_inc,o);
            discriminant = do^2-norm(o)^2+R^2;
            if discriminant >= 0
                s = -do - sqrt(discriminant);
                x_n = o + d_inc*s;

                n_vec = x_n;
                d(n,:) = d_inc-2*dot(d_inc,n_vec)*n_vec; % reflected ray

                O(n,:,2) = x_n;

        %         normals(n,:) = n_vec;
        %         O(n,:,3) = O(n,:,2) + 1.1*r;
                objectHit(n) = true;
        %         if plotRays
        %             hold on
        %             x = O(n,1,1:2);
        %             y = O(n,2,1:2);
        %             z = O(n,3,1:2);
        %             plot3(x(:),y(:),z(:),'black')
        %             x = [task.varCol{1}.o(n,1); reshape(O(n,1,2:end),size(O,3)-1,1)];
        %             y = [task.varCol{1}.o(n,2); reshape(O(n,2,2:end),size(O,3)-1,1)];
        %             z = [task.varCol{1}.o(n,3); reshape(O(n,3,2:end),size(O,3)-1,1)];
        %             plot3(x(:),y(:),z(:),'black')
        %         end
            end
        end
    case 'M3'
        Eps = 1e6*eps;
        L = task.varCol{1}.L;
        R1 = task.varCol{1}.R1;
        R2 = task.varCol{1}.R2;
        x0 = L*R2/(R1-R2);
        mu = (R2/x0)^2;
        o = O(:,:,1);
        d_INC = repmat(d_inc,N,1);
        do = dot(d_INC,o,2);
        s = Inf(N,3);
        discriminant = do.^2-norm2(o).^2+R2^2;
        stemp = -do - sqrt(discriminant);
        x_n = o + stemp*d_inc;
        indices = and(discriminant >= 0, x_n(:,1) >= -Eps);
        s(indices,1) = stemp(indices);
        a = 1;
        b = 2*d_inc(1)*L+2*do;
        c = norm2(o).^2 + 2*o(:,1)*L + L^2-R1^2;

        discriminant = b.^2-4*a*c;
        stemp = (-b-sqrt(discriminant))/(2*a);
        x_n = o + stemp*d_inc;
        indices = and(discriminant >= 0, x_n(:,1)+L <= Eps);
        s(indices,2) = stemp(indices);
        a = d_inc(2)^2+d_inc(3)^2-mu*d_inc(1)^2;
        b = 2*o(:,3)*d_inc(3)+2*o(:,2)*d_inc(2)-2*mu*o(:,1)*d_inc(1)+2*d_inc(1)*x0*mu;
        c = o(:,2).^2+o(:,3).^2 - mu*o(:,1).^2 + 2*mu*o(:,1)*x0 - mu*x0^2;

        discriminant = b.^2-4*a*c;
        stemp = (-b-sqrt(discriminant))/(2*a);
        x_n = o + stemp*d_inc;
        indices = and(discriminant >= 0, and(x_n(:,1) < 0, x_n(:,1) > -L));
        s(indices,3) = stemp(indices);
        objectHit = any(~isinf(s),2);

        [mins,I] = min(s,[],2);
        x_n = o + mins*d_inc;
        n_vec = zeros(size(o));
        indices = I == 1;
        n_vec(indices,:) = x_n(indices,:)./repmat(norm2(x_n(indices,:)),1,3);
        indices = I == 2;
        xtemp = x_n(indices,:);
        xtemp(:,1) = xtemp(:,1) + L;
        n_vec(indices,:) = xtemp./repmat(norm2(xtemp),1,3);
        indices = I == 3;
        theta = atan2(x_n(:,2),x_n(:,3));
        n_vectemp = [ones(N,1),sin(theta)/sqrt(mu),cos(theta)/sqrt(mu)];
        n_vec(indices,:) = n_vectemp(indices,:)./repmat(norm2(n_vectemp(indices,:)),1,3);

        dn = dot(d_INC,n_vec,2);
        d = d_INC-2*dn(:,[1,1,1]).*n_vec; % reflected ray

        O(:,:,2) = x_n;
end
% if plotRays
%     axis equal
%     [X,Y,Z] = sphere(1000);
% %     surf(R2*X,R2*Y,R2*Z, 'FaceColor', 1.5*[44 77 32]/255,'LineStyle','none')
% %     axis off
%     camlight
%     keyboard
% end

temp = zeros(N,1);
temp(objectHit) = 1:sum(objectHit);

beams = task.varCol{1}.beams;
beams = reshape(temp(beams(:)),size(beams,1),size(beams,2));
beams(any(~beams,2),:) = [];

task.varCol{1}.O = O(objectHit,:,:);
task.varCol{1}.o = task.varCol{1}.o(objectHit,:);
task.varCol{1}.d = d(objectHit,:);
dofs = size(task.varCol{1}.O,1);
task.dofs = dofs;
task.varCol{1}.beams = beams;
% task.varCol{1}.A = p_inc(O(objectHit,:,end-1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% XX = task.varCol{1}.o;
% XX2 = task.varCol{1}.O(:,1:3,1);
% d_vec = task.d_vec;
% XX_m = orthogonalTransform(XX, d_vec);
% XX2_m = orthogonalTransform(XX2, d_vec);
% YY_m = XX_m(:,2);
% XX_m = XX_m(:,1);
% YY2_m = XX2_m(:,2);
% XX2_m = XX2_m(:,1);
%    
% oXX_m = task.varCol{1}.oXX_m;
% oYY_m = task.varCol{1}.oYY_m;
% X_m = task.varCol{1}.X_m;
% convexHull = task.varCol{1}.convexHull;
% close all
% figure(1)
% hold on
% axis equal
% t = linspace(0,2*pi,10000);
% plot(cos(t),sin(t))
% plot(oXX_m,oYY_m,'*','color','cyan')
% plot(XX_m,YY_m,'*','color','blue')
% plot(XX2_m,YY2_m,'*','color','green')
% plot(X_m(convexHull,1),X_m(convexHull,2),'*-','color','red')
% % for i = 1:size(XX2_m,1)
% %     text(XX2_m(i),YY2_m(i),num2str(i))
% % end
% % for i = 1:size(beams,1)
% %     indices = beams(i,:);
% %     h = plot(XX2_m(indices([2,3,4,5,6,7,2])),YY2_m(indices([2,3,4,5,6,7,2])),'r');
% %     pause(0.5)
% %     set(h,'Visible','off')
% % end     
% keyboard
