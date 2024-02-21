clear all


% nurbs = read_g2('NURBSgeometries/g2files/S1patched.g2');
% nurbs = getEllipsoidData('parm',2,'t',0.1);
load('test.mat')
% nurbs = getEllipsoidData('parm',2);
% nurbs = getPrismData('d_p',2);

absError = 0;
for patch = 1:numel(nurbs)
    d = nurbs{patch}.d;
    d_p = nurbs{patch}.d_p;
    knots = nurbs{patch}.knots;
    degree = nurbs{patch}.degree;
    [Q, W] = gaussTensorQuad(10+max(degree));
    number = nurbs{patch}.number;
    for j = 1:d_p
        uniqueXi = unique(knots{j});
        i_uniqueXi = floor(1 + rand()*(numel(uniqueXi)-1));
        Xi_e = [uniqueXi(i_uniqueXi),uniqueXi(i_uniqueXi+1)];
        Xi_e(2) = Xi_e(1) + rand()*(Xi_e(2)-Xi_e(1));
        J_2 = (Xi_e(:,2)-Xi_e(:,1))/2;
        xi = [Xi_e(1); parent2ParametricSpace(Xi_e, Q); Xi_e(2)];
        I = zeros(size(xi,1),d_p);
        I(:,j) = findKnotSpan(number(j), degree(j), xi, knots{j});
        XI = zeros(size(xi,1),d_p);
        XI(:,j) = xi;
        for i = 1:d_p
            if i ~= j
                XI(:,i) = rand();
                I(:,i) = findKnotSpan(number(i), degree(i), XI(:,i), knots{i});
            end
        end
        weights = reshape(nurbs{patch}.coeffs(d+1,:),1,[]);
        
        order = degree+1;
        n_en = prod(order);
        npts = size(XI,1);
        
        w_i = ones([npts,order]);
        shift = 1;
        for i = 1:d_p
            temp_i = ones(1,d_p);
            temp_i(i) = degree(i)+1;
            w_i = w_i + (reshape(I(:,i)-degree(i)+(0:degree(i)),[npts,temp_i])-1)*shift;
            shift = shift*number(i);
        end
        w = weights(w_i);
        RdR = NURBSbasis(I,XI,degree,knots,w,2,1); % Compute all second derivatives
    
        switch d_p
            case 2
                R = RdR{1};
                dRdxi = RdR{2}(:,:,1);
                dRdeta = RdR{3}(:,:,1);
                d2Rdxi2 = RdR{2}(:,:,2);
                d2Rdeta2 = RdR{3}(:,:,2);
                d2Rdxieta = RdR{4};

                absError = absError + sum((R(end,:) - R(1,:) - sum(RdR{j+1}(2:end-1,:,1).*W*J_2,1)).^2)./sum(R.^2,'all');
                absError = absError + sum((RdR{j+1}(end,:,1) - RdR{j+1}(1,:,1) - sum(RdR{j+1}(2:end-1,:,2).*W*J_2,1)).^2)./sum(RdR{j+1}(:,:,1).^2,'all');
                switch j
                    case 1
                        absError = absError + sum((dRdeta(end,:) - dRdeta(1,:) - sum(d2Rdxieta(2:end-1,:).*W*J_2,1)).^2)./sum(dRdeta.^2,'all');
                    case 2
                        absError = absError + sum((dRdxi(end,:) - dRdxi(1,:) - sum(d2Rdxieta(2:end-1,:).*W*J_2,1)).^2)./sum(dRdxi.^2,'all');
                end
            case 3
                R = RdR{1};
                dRdxi = RdR{2}(:,:,1);
                dRdeta = RdR{3}(:,:,1);
                dRdzeta = RdR{4}(:,:,1);
                d2Rdetazeta = RdR{5};
                d2Rdxizeta = RdR{6};
                d2Rdxieta = RdR{7};

                absError = absError + sum((R(end,:) - R(1,:) - sum(RdR{j+1}(2:end-1,:,1).*W*J_2,1)).^2)./sum(R.^2,'all');
                absError = absError + sum((RdR{j+1}(end,:,1) - RdR{j+1}(1,:,1) - sum(RdR{j+1}(2:end-1,:,2).*W*J_2,1)).^2)./sum(RdR{j+1}(:,:,1).^2,'all');
                switch j
                    case 1
                        absError = absError + sum((dRdzeta(end,:) - dRdzeta(1,:) - sum(d2Rdxizeta(2:end-1,:).*W*J_2,1)).^2)./sum(dRdzeta.^2,'all');
                        absError = absError + sum((dRdeta(end,:) - dRdeta(1,:) - sum(d2Rdxieta(2:end-1,:).*W*J_2,1)).^2)./sum(dRdeta.^2,'all');
                    case 2
                        absError = absError + sum((dRdzeta(end,:) - dRdzeta(1,:) - sum(d2Rdetazeta(2:end-1,:).*W*J_2,1)).^2)./sum(dRdzeta.^2,'all');
                        absError = absError + sum((dRdxi(end,:) - dRdxi(1,:) - sum(d2Rdxieta(2:end-1,:).*W*J_2,1)).^2)./sum(dRdxi.^2,'all');
                    case 3
                        absError = absError + sum((dRdeta(end,:) - dRdeta(1,:) - sum(d2Rdetazeta(2:end-1,:).*W*J_2,1)).^2)./sum(dRdeta(1,:).^2,'all');
                        absError = absError + sum((dRdxi(end,:) - dRdxi(1,:) - sum(d2Rdxizeta(2:end-1,:).*W*J_2,1)).^2)./sum(dRdxi.^2,'all');
                end
                
        end
    end
end
sqrt(absError)