function task = computeROMsolution(task,printLog)

noDofs = size(task.V,1);
omega_P = task.rom.omega;
omega = task.misc.omega;
basisROM = task.rom.basisROM;
task.misc.printLog = false;

if ~task.rom.adaptiveROM
    noVecs = task.rom.noVecs;
    P = numel(omega_P);
    omega_start = omega_P(1);
    omega_end = omega_P(end);
    if printLog
        fprintf('\nBasis:%s, noVecs = %d', basisROM, noVecs)
    end
    if printLog
        fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Computing basis for ROM ... ')
    end
    t_startROM = tic;
    switch basisROM
        case 'DGP'
            task = buildDGP(task);
            task.V = task.P_rightinv*task.V; % Scale back to match scaling for task.U
            task = rmfields(task,{'A','FF','Pinv','A0','A1','A2','A4'});
        case 'Hermite'
    %         mp.Digits(400);
    %         Y = getInterpolatingHermite(mp(k_P.'),mp(omega_ROM),noVecs);
    %         Y = double(Y);
            Y = getInterpolatingHermite(omega_P.',omega,noVecs);
        case 'Pade'
            p = cell(P,1);
            q = cell(P,1);
            useHP = 0;
            if useHP
                mp.Digits(400);
            end
            for i = 1:P
                n = ceil(noVecs/2)-1;
                m = floor(noVecs/2);
                shift = (i-1)*noVecs;
                if useHP
                    f = @(x,n) mp(task.V(:,shift+n+1));
                    [p{i},q{i}] = pade(f,mp(omega_P(i)),n,m); 
                else
                    f = @(x,n) task.V(:,shift+n+1);
                    [p{i},q{i}] = pade(f,omega_P(i),n,m); 
                end
            end
        case 'Splines'
            Vsize = P*noVecs;
    %         p_ROM = noVecs;
            p_ROM = P*noVecs-1; % optimal conditioning and precision
            GLL = GLLpoints(Vsize-(p_ROM+1)+2);     
            Xi = [zeros(1,p_ROM+1),parent2ParametricSpace([0,1],GLL(2:end-1)),ones(1,p_ROM+1)];
            task.V = zeros(Vsize);
            b = zeros(Vsize,noDofs);
    
            dp = zeros(noVecs,Vsize,P);
            for i = 1:P
                xi = (omega_P(i)-omega_start)/(omega_end - omega_start);
                i1 = findKnotSpan(Vsize, p_ROM, xi, Xi);
                ders = Bspline_basisDers3(i1, xi, p_ROM, Xi, noVecs-1);
                for j = 2:noVecs
                    ders(j,:) = ders(j,:)/(omega_end - omega_start).^(j-1);
                end
                dp(:,i1-p_ROM:i1,i) = ders;
            end
            counter = 1;
            for i = 1:noVecs
                for j = 1:P
                    shift = (j-1)*noVecs;
                    task.V(counter,:) = dp(i,:,j);
                    b(counter,:) = task.V(:,shift+i).';
                    counter = counter + 1;
                end
            end
            cond(task.V)
            a = task.V\b;
        case 'Bernstein'
            Vsize = P*noVecs;
            p_ROM = P*noVecs-1;
            useHP = 0;
            if useHP
                mp.Digits(400);
                p_ROM = mp(p_ROM);
            end
            omega_P = convert(omega_P,class(p_ROM));
            omega_end = omega_P(end);
            omega_start = omega_P(1);
            tic
            B = bernsteinBasis((omega_P-omega_start)/(omega_end - omega_start),p_ROM,noVecs-1,1/(omega_end - omega_start));
            toc
            task.V = zeros(Vsize,class(B));
            b = zeros(Vsize,noDofs);
            counter = 1;
            for i = 1:noVecs
                for j = 1:P
                    shift = (j-1)*noVecs;
                    task.V(counter,:) = B(j,:,i);
                    b(counter,:) = task.V(:,shift+i).';
                    counter = counter + 1;
                end
            end
    %             invV = inv(V);
    %         cond(V)
    %         a = double(invV*mp(b(:,1)));
    %         a = double(invV*mp(b));
            tic
    %                 a = double(V\mp(b));
            a = double(task.V\b);
            toc
    %         a = double(invV)*b;
    
        case 'Lagrange'
            Vsize = P*noVecs;
            P = noTasks;
            n = Vsize;
            temp = cos(pi/(2*n));
            ak = ((omega_start+omega_end)*temp-(omega_end-omega_start))/(2*temp);
            bk = ((omega_start+omega_end)*temp+(omega_end-omega_start))/(2*temp);
            j = 1:n;
            k_arrLagr = 1/2*(ak+bk)+1/2*(bk-ak)*cos((2*n-2*j+1)*pi/(2*n));
    
            dp = zeros(P,noVecs,n);
            for i = 1:n
                dp(:,1,i) = lagrangePolynomials(omega_P.',i,n,k_arrLagr);
                dp(:,2:noVecs,i) = lagrangePolynomialsNthDeriv(omega_P.',i,n,k_arrLagr,noVecs-1);
            end
            task.V = zeros(Vsize);
            b = zeros(Vsize,noDofs);
            counter = 1;
            for i = 1:noVecs
                for j = 1:P
                    shift = (j-1)*noVecs;
                    temp = dp(j,i,:);
                    task.V(counter,:) = temp(:);
                    b(counter,:) = task.V(:,shift+i).';
                    counter = counter + 1;
                end
            end
            cond(task.V)
            a = task.V\b;
        case 'Fourier'
            Vsize = P*noVecs;
            n_arr = 1:Vsize-1;
            counter = 1;
            task.V = zeros(Vsize);
            b = zeros(Vsize,noDofs);
            task.V(1:P,1) = 1;
            for i = 1:noVecs
                for j = 1:P
                    shift = (j-1)*noVecs;
                    if mod(i,2)
                        task.V(counter,2:end) = (-1)^((i-1)/2)*(n_arr*pi/(omega_end-omega_start)).^(i-1).*cos(n_arr*pi*omega_P(j)/(omega_end-omega_start));
                    else
                        task.V(counter,2:end) = (-1)^(i/2)*(n_arr*pi/(omega_end-omega_start)).^(i-1).*sin(n_arr*pi*omega_P(j)/(omega_end-omega_start));
                    end
                    b(counter,:) = task.V(:,shift+i).';
                    counter = counter + 1;
                end
            end
            cond(task.V)
            a = task.V\b;
    end  
    if printLog
        fprintf('using %12f seconds.', toc(t_startROM))
    end
end
%     k_arr3 = linspace(k_start,k_end,100);
%     k_arr3 = sort(unique([k_P, k_arr3]));
%     nPts = numel(k_arr3);
% 
%     p = zeros(size(k_arr3));
%     switch basisROM
%         case 'Splines'
%             N = zeros(Vsize,nPts);
%             for i = 1:nPts
%                 xi = (k_arr3(i)-k_start)/(k_end - k_start);
%                 i1 = findKnotSpan(Vsize, p_ROM, xi, Xi);
%                 ders = Bspline_basisDers3(i1, xi, p_ROM, Xi, p_ROM);
%                 N(i1-p_ROM:i1,i) = ders(1,:);
%             end
%             p = a(:,1).'*N;
%         case 'Lagrange'
%             for i = 1:Vsize
%                 p = p + a(i,1)*lagrangePolynomials(k_arr3.',i,n,k_arrLagr).';
%             end
%         case 'Taylor'
%             for i = 1:Vsize
%                 p = p + a(i,1)*(k_arrLagr-k_m).^(i-1)./factorial(i-1);
%             end
%         case 'Fourier'
%             for i = 1:Vsize
%                 p = p + a(i,1)*cos((i-1)*pi*k_arr3/(k_end-k_start));
%             end
%         case 'Bernstein'
%             B = bernsteinBasis((k_arr3-k_start)/(k_end - k_start),p_ROM,0);
%             p = a(:,1).'*B.';
%     end
% k_arr = double(k_arr);
% e3Dss_options.omega = 1500*k_arr;
%     x = evaluateNURBS(fluid{1},[0,0,0]).';
% p_ref = analytic_(x,e3Dss_options);
% close all
% semilogy(k_arr,abs(p-p_ref)./abs(p_ref),'DisplayName',basisROM)
% % hold on
% % uiopen('C:\Users\Zetison\Dropbox\work\matlab\ROMlagrange.fig',1)
% legend show
% figure(42)
% plot(k_arr,real(p),k_arr,real(p_ref))

% p = [U_P{1}(1,1), U_P{2}(1,1), U_P{3}(1,1)];
% dp = [U_P{1}(1,2), U_P{2}(1,2), U_P{3}(1,2)];
% e3Dss_options.omega = 1500*double(k_P);
% p_ref = analytic_(x,e3Dss_options);
% abs(p-p_ref)./abs(p_ref)
% R = 1;
% [p_ref,dp_ref] = rigidSphereScattering(R,pi,double(k_P),R,1,100);
% abs(dp-dp_ref)./abs(dp_ref)
% 

%             omega_ROM = double(unique(sort([k_P,omega_ROM])));
omega = double(omega);
energyError = zeros(1,size(omega,2));
L2Error = zeros(1,size(omega,2));
H1Error = zeros(1,size(omega,2));
H1sError = zeros(1,size(omega,2));
surfaceError = zeros(1,size(omega,2));
p_h = zeros(1,size(omega,2));
for i_f = 1:numel(omega)
    omega_i = omega(i_f);
    task.misc.omega = omega_i;
    if printLog
        fprintf(['\n%-' num2str(task.misc.stringShift) 's'], ['Computing ROM solution (' num2str(i_f) '/' num2str(numel(omega)) ')... '])
    end
    t_startROM = tic;
    task = getAnalyticSolutions(task);
    switch basisROM
        case 'DGP'
            task = buildRHS(task,true);
            task = collectMatrices(task,false,true,false);
            FFm = task.V'*task.FF;  

            Am = task.A0_am + omega_i*task.A1_am + omega_i^2*task.A2_am + omega_i^4*task.A4_am;
            if strcmp(task.sol.preconditioner,'none')
                task.UU = task.V*(Am\FFm);
            else
                Pinv = spdiags(1./sqrt(diag(Am)),0,size(Am,1),size(Am,2));
                task.UU = task.V*(Pinv*((Pinv*Am*Pinv)\(Pinv*FFm)));
            end
        case 'Hermite'
            task.UU = zeros(noDofs,1);
            counter = 1;
            for i = 1:P
                shift = (i-1)*noVecs;
                for n = 1:noVecs
                    task.UU = task.UU + task.V(:,shift+n)*Y(counter,i_f);
                    counter = counter + 1;
                end
            end
        case 'Pade'
            if useHP
                task.UU = double(interPade(mp(omega_i),mp(omega_P),p,q));
            else
                task.UU = interPade(omega_i,omega_P,p,q);
            end
        case 'Bernstein'
            B = bernsteinBasis(double((omega_i-omega_start)/(omega_end - omega_start)),double(p_ROM),0);
            task.UU = (B*a).';
        case 'Taylor'
            task.UU = interTaylor(omega_i,omega_P,task.V,noVecs-1);
    end
    task = computeDerivedQuantities(task);
    task = postProcessSolution(task);
    if printLog
        fprintf('using %12f seconds.', toc(t_startROM))
    end

    if task.err.calculateSurfaceError || task.err.calculateVolumeError
        task = calculateErrors(task, 1, task.misc.stringShift);
        if task.err.calculateSurfaceError
            surfaceError(i_f) = task.results.surfaceError;
        end
        if task.err.calculateVolumeError
            L2Error(i_f) = task.results.L2Error;
            H1Error(i_f) = task.results.H1Error;
            H1sError(i_f) = task.results.H1sError;
            energyError(i_f) = task.results.energyError;
        end
    end

    if task.ffp.calculateFarFieldPattern
        task = calculateTS(task,printLog,task.misc.stringShift);
        p_h(i_f) = task.ffp.p_h;
    end
end
if task.ffp.calculateFarFieldPattern
    task.misc.omega = omega;
    task = getAnalyticSolutions(task);
    task = computeDerivedFFPquantities(task,p_h);
    fieldCell = {'p','p_Re','p_Im','abs_p','TS'};
    if task.analyticSolutionExist
        fieldCell = [fieldCell, 'p_ref', 'p_Re_ref', 'p_Im_ref','abs_p_ref','TS_ref','error_pAbs','error_p'];
    end
    for field = fieldCell
        switch basisROM
            case {'Taylor','Pade'}
                [task.results.(field{1}),temp_omega] = insertNaN(omega,task.results.(field{1}),omega_P);
            otherwise
                task.results.(field{1}) = task.results.(field{1});
                temp_omega = omega;
        end
    end
end
if task.err.calculateSurfaceError || task.err.calculateVolumeError
    switch basisROM
        case {'Taylor','Pade'}
            if task.err.calculateSurfaceError
                [surfaceError,temp_omega] = insertNaN(omega,surfaceError,omega_P);
            end
            if task.err.calculateVolumeError
                [energyError,temp_omega] = insertNaN(omega,energyError,omega_P);
                L2Error  = insertNaN(omega,L2Error,omega_P);
                H1Error  = insertNaN(omega,H1Error,omega_P);
                H1sError = insertNaN(omega,H1sError,omega_P);
            end
        otherwise
            temp_omega = omega;
    end
    if task.err.calculateVolumeError
        task.results.energyError = energyError;
        task.results.L2Error = L2Error;
        task.results.H1Error = H1Error;
        task.results.H1sError = H1sError;
    end
    if task.err.calculateSurfaceError
        task.results.surfaceError = surfaceError;
    end
end

task = rmfield(task,{'V'});
task.varCol = rmfields(task.varCol,getAddedFields());
if task.misc.clearGlobalMatrices
    task = rmfields(task,{'A','FF','Pinv','A0','A1','A2','A4','A0_am','A1_am','A2_am','A4_am'});
end
task.omega = temp_omega;
task.f = temp_omega/(2*pi);
task.varCol{1}.k = temp_omega/task.varCol{1}.c_f;

task.misc.printLog = true;

function [y,x] = insertNaN(x,y,a)
inter = (a(2:end)+a(1:end-1))/2;
for i = 1:numel(a)-1
    [x,I] = sort([x, inter(i)]);
    y = [y, NaN];
    y = y(I);
end

