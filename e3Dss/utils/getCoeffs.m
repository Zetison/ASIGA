function CC = getCoeffs(n, omega, a, b, k, options)

% Define the system size, I
if n == 0
    I = 4;
else
    I = 6;
end
ESBC = options.ESBC;
SHBC = options.SHBC;
SSBC = options.SSBC;
rho_f = options.rho_f;
P_inc = options.P_inc;
R_i = options.R_i;
R_o = options.R_o;
M = length(R_o);
r_s = options.r_s;

if ~(SHBC && M == 1)
    G = options.G;
end

if isa(R_o,'sym')
    PI = vpa('pi');
elseif isa(R_o,'mp')
    PI = mp('pi');
else
    PI = pi;
end

%% Calculate submatrices
H1 = zeros(I-2,I-2,length(omega),M,class(R_o));
H2 = zeros(2,length(omega),M,2,class(R_o));
H3 = zeros(2,length(omega),M,2,class(R_o));
H4 = zeros(I-2,length(omega),M,2,class(R_o));

H21 = zeros(length(omega),1,class(R_o));
D2 = zeros(length(omega),1,class(R_o));

if SHBC && M == 1 % as this case is independent of rho_f, we here avoid demanding its value
    rho_f_omega2 = omega.^2;
else
    rho_f_omega2 = rho_f(1)*omega.^2;
end
for m = 1:M+1
    if m == 1
        zeta = k(:,1)*R_o(1);
        H11 = -1./rho_f_omega2.*(n*hankel_s(n,zeta,1) - zeta.*hankel_s(n+1,zeta,1));
        if options.usePlaneWave
            if isa(P_inc,'function_handle') 
                D1 = P_inc(omega)*(2*n+1)*1i^n./rho_f_omega2.*dbessel_s(n,zeta,1,[],1);  % dbessel_s = zeta*dj_n
            else
                D1 = P_inc*(2*n+1)*1i^n./rho_f_omega2.*dbessel_s(n,zeta,1,[],1);  % dbessel_s = zeta*dj_n
            end
        elseif options.usePointChargeWave
            q = @(v) sqrt(R_o(1)^2 + 2*r_s*R_o(1)*v + r_s^2);
            F2 = zeros(size(k,1),1,class(R_o));
            
            parfor i = 1:size(k,1)
                Phi_k = @(v) exp(1i*k(i,1)*q(v))./(4*PI*q(v));
                integrand = @(v) 4*PI*r_s*(R_o(1) + r_s*v).*Phi_k(v)./q(v).^2.*(1i*k(i,1)*q(v) - 1).*legendre_(n,v);
                if isa(P_inc,'function_handle')
                    F2(i) = P_inc(omega(i))*(2*n+1)/2*integral(integrand,-1,1);
                else
                    F2(i) = P_inc*(2*n+1)/2*integral(integrand,-1,1);
                end
            end
            D1 = R_o(1)./rho_f_omega2.*F2;
        end
        if ~(SHBC && M == 1)
            H21 = R_o(1)^2/(2*G(1))*hankel_s(n,zeta,1);   
            if options.usePlaneWave    
                if isa(P_inc,'function_handle') 
                    D2 = -P_inc(omega)*R_o(1)^2/(2*G(1))*(2*n+1)*1i^n.*bessel_s(n,k(:,1)*R_o(1),1);
                else
                    D2 = -P_inc*R_o(1)^2/(2*G(1))*(2*n+1)*1i^n*bessel_s(n,k(:,1)*R_o(1),1);
                end
            elseif options.usePointChargeWave
                F1 = zeros(size(k,1),1,class(R_o));
                parfor i = 1:size(k,1)
                    Phi_k = @(v) exp(1i*k(i,1)*q(v))./(4*PI*q(v));
                    integrand = @(v) 4*PI*r_s*Phi_k(v).*legendre_(n,v);
                    if isa(P_inc,'function_handle')
                        F1(i) = P_inc(omega(i))*(2*n+1)/2*integral(integrand,-1,1);
                    else
                        F1(i) = P_inc*(2*n+1)/2*integral(integrand,-1,1);
                    end
                end
                D2 = -R_o(1)^2/(2*G(1))*F1;
            end
            if M == 1 && ESBC
                H1(:,:,:,1) = H1_solid_(n,a(:,1),b(:,1),R_o(1));
                H4(:,:,1,2) = H4_solid_(n,a(:,1),b(:,1),R_o(1)); % H4 at R_o
            else
                H1(:,:,:,1) = H1_(n,a(:,1),b(:,1),R_o(1),R_i(1));
                H4(:,:,1,2) = H4_(n,a(:,1),b(:,1),R_o(1)); % H4 at R_o
            end
            if ~(M == 1 && (SSBC || ESBC))
                H4(:,:,1,1) = H4_(n,a(:,1),b(:,1),R_i(1)); % H4 at R_i
            end
        end        
    elseif m < M+1
        H2(:,:,m,1) = H2_(n,G(m-1),k(:,m),R_i(m-1)); % H2 at R_i
        H3(:,:,m,1) = H3_(n,rho_f(m),k(:,m),omega,R_i(m-1)); % H3 at R_i
        H3(:,:,m,2) = H3_(n,rho_f(m),k(:,m),omega,R_o(m)); % H3 at R_o
        if ~(SHBC && M == m)
            H2(:,:,m,2) = H2_(n,G(m),k(:,m),R_o(m)); % H2 at R_o

            if M == m && ESBC
                H1(:,:,:,m) = H1_solid_(n,a(:,m),b(:,m),R_o(m));
                H4(:,:,m,2) = H4_solid_(n,a(:,m),b(:,m),R_o(m)); % H4 at R_o
            else
                H1(:,:,:,m) = H1_(n,a(:,m),b(:,m),R_o(m),R_i(m));
                H4(:,:,m,2) = H4_(n,a(:,m),b(:,m),R_o(m)); % H4 at R_o
            end
            if ~(M == m && (SSBC || ESBC))
                H4(:,:,m,1) = H4_(n,a(:,m),b(:,m),R_i(m)); % H4 at R_i
            end
        end
    else
        if ~(SSBC || ESBC || SHBC)
            zeta = k(:,M+1)*R_i(M);
            H2iMp1 = R_i(M)^2/(2*G(M))*bessel_s(n,zeta,1);
            H3iMp1 = -1./(rho_f(M+1)*omega.^2).*dbessel_s(n,zeta,1,[],1); % dbessel_s = zeta*dj_n
        else
            H2iMp1 = NaN(length(omega),1);
            H3iMp1 = NaN(length(omega),1);
        end        
    end
end

%% Calculate coefficients CC for each frequency
% Allocate memory for coefficients, right hand side and global matrix
if ESBC
    if n == 0
        systemSize = I*M-2; 
    else
        systemSize = I*M-3;   
    end
elseif SHBC
    systemSize = I*M-(I-1); 
elseif SSBC
    systemSize = I*M-1;   
else
    systemSize = I*M;
end
CC = zeros(length(omega),6*M,class(R_o));

% Loop over frequencies
for i = 1:length(omega)
% parfor i = 1:length(omega)
%     H = sparse([], [], [], systemSize,systemSize, M*((I-2)^2+2*(I-2)+8)); % global matrix
    H = zeros(systemSize,class(R_o)); % global matrix
    D = zeros(systemSize,1,class(R_o)); % righ hand side
    
    for m = 1:M+1
        if m == 1            
            H(1,1) = H11(i);
            D(1) = D1(i);
            if ~(SHBC && M == 1)
                H(2,1) = H21(i);
                D(2) = D2(i);
                if ESBC && M == 1
                    H4_temp = H4(:,i,1,2); % H4 at R_o
                    H1_temp = H1(:,:,i,1);

                    if n == 0
                        H(I*(M-1)+1,2) = H4_temp(1);
                        H(I*(M-1)+2,2) = H1_temp(1,1);
                    else
                        H(I*(M-1)+1,2:I-3) = H4_temp([1,3]);
                        H(I*(M-1)+2:I*M-3,2:I-3) = H1_temp([1,2],[1,3]);
                    end
                else
                    H(1:1,2:I-1) = H4(:,i,1,2); % H4 at R_o
                    H(2:I-1,2:I-1) = H1(:,:,i,1);
                end
                if ~(1 == M && (SSBC || ESBC))
                    H(I,2:I-1) = H4(:,i,1,1);  % H4 at R_i
                end  
            end

        elseif m < M+1
            H(I*(m-1)-1,I*(m-1):I*(m-1)+1) = H2(:,i,m,1);  % H2 at R_i(m-1)
            H(I*(m-1),I*(m-1):I*(m-1)+1) = H3(:,i,m,1);  % H3 at R_i(m-1)

            H(I*(m-1)+1,I*(m-1):I*(m-1)+1) = H3(:,i,m,2);  % H3 at R_o(m)
            if ~(SHBC && M == m)
                H(I*(m-1)+2,I*(m-1):I*(m-1)+1) = H2(:,i,m,2);  % H2 at R_o(m)
                if ESBC && M == m
                    H4_temp = H4(:,i,m,2);
                    H1_temp = H1(:,:,i,m);

                    if n == 0
                        H(I*(M-1)+1,I*(m-1)+2) = H4_temp(1);
                        H(I*(M-1)+2,I*(m-1)+2) = H1_temp(1,1);
                    else
                        H(I*(M-1)+1,I*(m-1)+2:I*m-3) = H4_temp([1,3]);
                        H(I*(M-1)+2:I*M-3,I*(m-1)+2:I*m-3) = H1_temp([1,2],[1,3]);
                    end
                else
                    H(I*(m-1)+1:I*(m-1)+1,I*(m-1)+2:I*m-1) = H4(:,i,m,2); % K4 at R_o
                    H(I*(m-1)+2:I*m-1,I*(m-1)+2:I*m-1) = H1(:,:,i,m);
                end
                if ~(m == M && (SSBC || ESBC))
                    H(I*m,I*(m-1)+2:I*m-1) = H4(:,i,m,1);  % K4 at R_i
                end 
            end
        else
            if ~(SSBC || ESBC || SHBC)
                H(I*M-1,I*M) = H2iMp1(i); 
                H(I*M,I*M) = H3iMp1(i);     
            end
        end
    end
%     keyboard
    Pinv = diag(1./max(abs(H)));
    if any(isinf(Pinv(:)))
        error('K was singular')
    end
    if n == 0
        indices = setdiff(1:6*M, sort([(4:6:6*M) (5:6:6*M)]));
    else
        indices = 1:6*M;
    end
    if ESBC
        if n == 0
            indices = indices(1:end-2);    
        else
            indices = indices([1:end-4,6*M-2]);    
        end
    elseif SHBC
        indices = indices(1:end-(I-1));   
    elseif SSBC
        indices = indices(1:end-1); 
    end
    % Fenders routine for reducing the complex system of linear equations
    % to a real system of linear equations is not implemented
    CC_temp = zeros(1,6*M,class(R_o));
    CC_temp(indices) = diag(Pinv).*((H*Pinv)\D);
    CC(i,:) = CC_temp;
    % Uncomment the following to get the spy matrix in the paper
%     if n == 300
%         keyboard
%         fileName = 'results/spy_H';
%         figure(1)
%         spy2(H)
%         cond(full(H))
%         extraAxisOptions = {...
%             'axis on top=true', ...
%             'at={(0,0)}', ...
%             'xtick={1,2,...,18}', ...
%             'ytick={1,2,...,18}', ...
%             'xlabel=$j$', ...
%             'ylabel=$i$', ...
%             'colorbar style={ylabel={$H_{ij,300}$}, ytick={-300,-200,...,200}, yticklabels={$10^{-300}$, $10^{-200}$, $10^{-100}$, $10^{0}$, $10^{100}$, $10^{200}$}}'};
%         
%         matlab2tikz([fileName '_1.tex'], 'height', '3.2094in', ...
%             'extraAxisOptions', extraAxisOptions, 'relativeDataPath', '../../../matlab/otherFunctions/e3Dss/results/') % , 'imagesAsPng', false
%         figure(2)
%         spy2(K,1)
%         cond(full(K))
%         matlab2tikz([fileName '_1_bw.tex'], 'height', '3.2094in', ...
%             'extraAxisOptions', extraAxisOptions, 'relativeDataPath', '../../../matlab/otherFunctions/e3Dss/results/')
%         figure(3)
%         spy2(K*Pinv)
%         cond(full(K*Pinv))
%         extraAxisOptions{7} = 'colorbar style={ylabel={$\tilde{H}_{ij,300}$}, ytick={-10,-8,...,0}, yticklabels={$10^{-10}$, $10^{-8}$, $10^{-6}$, $10^{-4}$, $10^{-2}$, $10^0$}}';
%         matlab2tikz([fileName '_2.tex'], 'height', '3.2094in', ...
%             'extraAxisOptions', extraAxisOptions, 'relativeDataPath', '../../../matlab/otherFunctions/e3Dss/results/')
%         figure(4)
%         spy2(K*Pinv,1)
%         cond(full(K*Pinv))
%         matlab2tikz([fileName '_2_bw.tex'], 'height', '3.2094in', ...
%             'extraAxisOptions', extraAxisOptions, 'relativeDataPath', '../../../matlab/otherFunctions/e3Dss/results/')
%     end
end

function H = H1_solid_(n,a,b,R_o)
xi = a*R_o;
eta = b*R_o;
Z{1,1} = bessel_s(n,xi,1);
Z{1,2} = bessel_s(n+1,xi,1);
if n == 0
    H = zeros(2,2,length(a),class(R_o));

    H(1,1,:) = S_(5, 1, n, xi, eta, Z);
else
    H = zeros(4,4,length(a),class(R_o));

    H(1,1,:) = S_(5, 1, n, xi, eta, Z);
    H(2,1,:) = S_(7, 1, n, xi, eta, Z);

    Z{1,1} = bessel_s(n,eta,1);
    Z{1,2} = bessel_s(n+1,eta,1);

    H(1,3,:) = T_(5, 1, n, eta, Z);
    H(2,3,:) = T_(7, 1, n, eta, Z);
end

function H = H1_(n,a,b,R_o,R_i)

xi = a*R_o;
eta = b*R_o;
Z{2,1} = bessel_s(n,xi,2);
Z{2,2} = bessel_s(n+1,xi,2);
Z{1,1} = bessel_s(n,xi,1);
Z{1,2} = bessel_s(n+1,xi,1);
if n == 0
    H = zeros(2,2,length(a),class(R_o));
    
    H(1,1,:) = S_(5, 1, n, xi, eta, Z);
    H(1,2,:) = S_(5, 2, n, xi, eta, Z);

    xi = a*R_i;
    eta = b*R_i;
    Z{1,1} = bessel_s(n,xi,1);
    Z{2,1} = bessel_s(n,xi,2);
    Z{1,2} = bessel_s(n+1,xi,1);
    Z{2,2} = bessel_s(n+1,xi,2);

    H(2,1,:) = S_(5, 1, n, xi, eta, Z);
    H(2,2,:) = S_(5, 2, n, xi, eta, Z);
else
    H = zeros(4,4,length(a),class(R_o));

    H(1,1,:) = S_(5, 1, n, xi, eta, Z);
    H(1,2,:) = S_(5, 2, n, xi, eta, Z);
    H(2,1,:) = S_(7, 1, n, xi, eta, Z);
    H(2,2,:) = S_(7, 2, n, xi, eta, Z);

    Z{1,1} = bessel_s(n,eta,1);
    Z{2,1} = bessel_s(n,eta,2);
    Z{1,2} = bessel_s(n+1,eta,1);
    Z{2,2} = bessel_s(n+1,eta,2);

    H(1,3,:) = T_(5, 1, n, eta, Z);
    H(1,4,:) = T_(5, 2, n, eta, Z);
    H(2,3,:) = T_(7, 1, n, eta, Z);
    H(2,4,:) = T_(7, 2, n, eta, Z);
    
    xi = a*R_i;
    eta = b*R_i;
    Z{1,1} = bessel_s(n,xi,1);
    Z{2,1} = bessel_s(n,xi,2);
    Z{1,2} = bessel_s(n+1,xi,1);
    Z{2,2} = bessel_s(n+1,xi,2);

    H(3,1,:) = S_(7, 1, n, xi, eta, Z);
    H(3,2,:) = S_(7, 2, n, xi, eta, Z);
    H(4,1,:) = S_(5, 1, n, xi, eta, Z);
    H(4,2,:) = S_(5, 2, n, xi, eta, Z);

    Z{1,1} = bessel_s(n,eta,1);
    Z{2,1} = bessel_s(n,eta,2);
    Z{1,2} = bessel_s(n+1,eta,1);
    Z{2,2} = bessel_s(n+1,eta,2);

    H(3,3,:) = T_(7, 1, n, eta, Z);
    H(3,4,:) = T_(7, 2, n, eta, Z);
    H(4,3,:) = T_(5, 1, n, eta, Z);
    H(4,4,:) = T_(5, 2, n, eta, Z);
end

function H = H2_(n,mu,k,R)

H = zeros(2,length(k),class(R));

H(1,:) = R^2/(2*mu)*bessel_s(n,k*R,1);
H(2,:) = R^2/(2*mu)*bessel_s(n,k*R,2);

function H = H3_(n,rho_f,k,omega,R)

H = zeros(2,length(k),class(R));

zeta = k*R;
H(1,:) = -1./(rho_f*omega.^2).*dbessel_s(n,zeta,1,[],1); % dbessel_s = zeta*dj_n
H(2,:) = -1./(rho_f*omega.^2).*dbessel_s(n,zeta,2,[],1); % dbessel_s = zeta*dy_n

function H = H4_solid_(n,a,b,R)

xi = a*R;
eta = b*R;
Z = cell(2,2);
Z{1,1} = bessel_s(n,xi,1);
Z{1,2} = bessel_s(n+1,xi,1);
if n == 0
    H = zeros(2,length(a),class(R));
    
    H(1,:) = S_(1, 1, n, xi, eta, Z);
else
    H = zeros(4,length(a),class(R));
    
    H(1,:) = S_(1, 1, n, xi, eta, Z);
    
    Z{1,1} = bessel_s(n,eta,1);
    Z{1,2} = bessel_s(n+1,eta,1);
    H(3,:) = T_(1, 1, n, eta, Z);
end


function H = H4_(n,a,b,R)

xi = a*R;
eta = b*R;
Z = cell(2,2);
Z{1,1} = bessel_s(n,xi,1);
Z{1,2} = bessel_s(n+1,xi,1);
Z{2,1} = bessel_s(n,xi,2);
Z{2,2} = bessel_s(n+1,xi,2);

if n == 0
    H = zeros(2,length(a),class(R));
    
    H(1,:) = S_(1, 1, n, xi, eta, Z);
    H(2,:) = S_(1, 2, n, xi, eta, Z);
else
    H = zeros(4,length(a),class(R));
    
    H(1,:) = S_(1, 1, n, xi, eta, Z);
    H(2,:) = S_(1, 2, n, xi, eta, Z);
    
    Z{1,1} = bessel_s(n,eta,1);
    Z{1,2} = bessel_s(n+1,eta,1);
    H(3,:) = T_(1, 1, n, eta, Z);
    Z{2,1} = bessel_s(n,eta,2);
    Z{2,2} = bessel_s(n+1,eta,2);
    H(4,:) = T_(1, 2, n, eta, Z);
end








