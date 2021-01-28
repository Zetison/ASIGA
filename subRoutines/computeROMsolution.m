function tasks = computeROMsolution(tasks,i_task,basisROMcell,omega_ROM,noVecsArr)
% basisROM = 'Taylor';
% basisROM = 'Pade';
% basisROM = 'Lagrange';
% basisROM = 'Splines';
% basisROM = 'Fourier';
% basisROM = 'Bernstein';
runTasksInParallel = false;

U_P = tasks(i_task).task.varCol{1}.U_sweep;
U_P2 = tasks(i_task).task.varCol{1}.U_sweep2;
noDofs = size(U_P{1},1);
noDofs2 = size(U_P2{1},1);
A0 = tasks(i_task).task.varCol{1}.A0;
A1 = tasks(i_task).task.varCol{1}.A1;
A2 = tasks(i_task).task.varCol{1}.A2;
tasks(i_task).task.varCol = rmfields(tasks(i_task).task.varCol,{'A0','A1','A2','U_sweep','U_sweep2','U'});
task = tasks(i_task).task;
omega_P = task.omega;
P = numel(omega_P);
omega_start = omega_P(1);
omega_end = omega_P(end);
stringShift = 40;
varCol = task.varCol;
tasks(i_task).task.varCol = extractVarColFields(task,varCol);
noDomains = numel(varCol);

% noVec = size(U_P{1},2);
for i_b = 1:numel(basisROMcell)
    basisROM = basisROMcell{i_b};
    for taskROM = 1:numel(noVecsArr)
        noVecs = noVecsArr(taskROM);

        fprintf(['\n%-' num2str(stringShift) 's'], 'Computing basis for ROM ... ')
        t_startROM = tic;
        switch basisROM
            case 'DGP'
                V = zeros(noDofs2,P*noVecs);
                counter = 1;
                for i = 1:noVecs
                    for j = 1:P
                        V(:,counter) = U_P2{j}(:,i);
                        counter = counter + 1;
                    end
                end
                allDofsToRemove = varCol{1}.allDofsToRemove;
                V(allDofsToRemove,:) = [];
                U = zeros(size(V));
                U(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
                for i = 2:size(V,2)
                    U(:,i) = V(:,i);
                    for j = 1:i-1
                        U(:,i) = U(:,i) - ( U(:,j)'*U(:,i) )/( U(:,j)'*U(:,j) )*U(:,j);
                    end
                    U(:,i) = U(:,i)/sqrt(U(:,i)'*U(:,i));
                end
                V = U;
                A0_am = V'*A0*V;
                A1_am = V'*A1*V;
                A2_am = V'*A2*V;
            case 'Hermite'
%                 mp.Digits(400);
%                 Y = getInterpolatingHermite(mp(k_P.'),mp(omega_ROM),noVecs);
%                 Y = double(Y);
                Y = getInterpolatingHermite(omega_P.',omega_ROM,noVecs);
            case 'Pade'
                p = cell(P,1);
                q = cell(P,1);
                useHP = 0;
                if useHP
                    mp.Digits(400);
%                     mp.Digits(1000);
                end
                for i = 1:P
                    n = ceil(noVecs/2)-1;
                    m = floor(noVecs/2);
                    if useHP
                        f = @(x,n) mp(U_P{i}(:,n+1));
                        [p{i},q{i}] = pade(f,mp(omega_P(i)),n,m); 
                    else
                        f = @(x,n) U_P{i}(:,n+1);
                        [p{i},q{i}] = pade(f,omega_P(i),n,m); 
                    end
                end
                
            case 'Splines'
                Vsize = P*noVecs;
        %         p_ROM = noVecs;
                p_ROM = P*noVecs-1; % optimal conditioning and precision
                GLL = GLLpoints(Vsize-(p_ROM+1)+2);     
                Xi = [zeros(1,p_ROM+1),parent2ParametricSpace([0,1],GLL(2:end-1)),ones(1,p_ROM+1)];
                V = zeros(Vsize);
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
                        V(counter,:) = dp(i,:,j);
                        b(counter,:) = U_P{j}(:,i).';
                        counter = counter + 1;
                    end
                end
                cond(V)
                a = V\b;
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
                V = zeros(Vsize,class(B));
                b = zeros(Vsize,noDofs);
                counter = 1;
                for i = 1:noVecs
                    for j = 1:P
                        V(counter,:) = B(j,:,i);
                        b(counter,:) = U_P{j}(:,i).';
                        counter = counter + 1;
                    end
                end
    %             invV = inv(V);
        %         cond(V)
        %         a = double(invV*mp(b(:,1)));
        %         a = double(invV*mp(b));
                tic
%                 a = double(V\mp(b));
                a = double(V\b);
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
                V = zeros(Vsize);
                b = zeros(Vsize,noDofs);
                counter = 1;
                for i = 1:noVecs
                    for j = 1:P
                        temp = dp(j,i,:);
                        V(counter,:) = temp(:);
                        b(counter,:) = U_P{j}(:,i).';
                        counter = counter + 1;
                    end
                end
                cond(V)
                a = V\b;
            case 'Fourier'
                Vsize = P*noVecs;
                n_arr = 1:Vsize-1;
                counter = 1;
                V = zeros(Vsize);
                b = zeros(Vsize,noDofs);
                V(1:P,1) = 1;
                for i = 1:noVecs
                    for j = 1:P
                        if mod(i,2)
                            V(counter,2:end) = (-1)^((i-1)/2)*(n_arr*pi/(omega_end-omega_start)).^(i-1).*cos(n_arr*pi*omega_P(j)/(omega_end-omega_start));
                        else
                            V(counter,2:end) = (-1)^(i/2)*(n_arr*pi/(omega_end-omega_start)).^(i-1).*sin(n_arr*pi*omega_P(j)/(omega_end-omega_start));
                        end
                        b(counter,:) = U_P{j}(:,i).';
                        counter = counter + 1;
                    end
                end
                cond(V)
                a = V\b;
        end  
        fprintf('using %12f seconds.', toc(t_startROM))
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

        calculateSurfaceError = 1;
        if calculateSurfaceError
%             omega_ROM = double(unique(sort([k_P,omega_ROM])));
            omega_ROM = double(omega_ROM);
            fprintf(['\n%-' num2str(stringShift) 's'], 'Computing ROM solution ... ')
            t_startROM = tic;
            switch basisROM
                case 'DGP'
                    U_fluid_oArr = zeros(noDofs2,numel(omega_ROM));
                    varCol{1}.omega = omega_ROM;
                    varCol{1}.k = omega_ROM/varCol{1}.c_f;
                    varCol{1}.noRHSs = numel(omega_ROM);
                    varCol = getAnalyticSolutions(varCol);
                    if noDomains > 1
                        varCol{2}.p_inc_ROM = varCol{1}.p_inc_ROM;
                        varCol{2}.dp_inc_ROM = varCol{1}.dp_inc_ROM;
                        varCol{2}.noRHSs = varCol{1}.noRHSs;
                    end
                    for i = 1:min(noDomains,2)
                        varCol{i}.useROM = false;
                        varCol{i} = applyNeumannCondition(varCol{i},i==2);
                        varCol{i}.useROM = true;
                    end
                    [varCol,FF] = collectMatrices(varCol,task);
                    FF = V'*FF;  
                    freeDofs = setdiff(1:noDofs2,allDofsToRemove);
                    
                    for i_f = 1:numel(omega_ROM)
                        omega = omega_ROM(i_f);
                        Am = A0_am + omega*A1_am + omega^2*A2_am;
                        U_fluid_oArr(freeDofs,i_f) = V*(Am\FF(:,i_f));
                        U_fluid_oArr(:,i_f) = addSolutionToRemovedNodes_new(U_fluid_oArr(:,i_f), varCol{1});
                    end
                    U_fluid_oArr(noDofs+1:end,:) = [];
                case 'Hermite'
                    U_fluid_oArr = zeros(noDofs,numel(omega_ROM));
                    counter = 1;
                    for i = 1:P
                        for n = 1:noVecs
                            U_fluid_oArr = U_fluid_oArr + U_P{i}(:,n)*Y(counter,:);
                            counter = counter + 1;
                        end
                    end
                case 'Pade'
                    if useHP
                        U_fluid_oArr = double(interPade(mp(omega_ROM),mp(omega_P),p,q));
                    else
                        U_fluid_oArr = interPade(omega_ROM,omega_P,p,q);
                    end
                case 'Bernstein'
                    B = bernsteinBasis(double((omega_ROM-omega_start)/(omega_end - omega_start)),double(p_ROM),0);
                    U_fluid_oArr = (B*a).';
                case 'Taylor'
                    U_fluid_oArr = interTaylor(omega_ROM,omega_P,U_P,noVecs-1);
            end
            fprintf('using %12f seconds.', toc(t_startROM))
            
            fprintf(['\n%-' num2str(stringShift) 's'], 'Computing errors for ROM sweeps ... ')
            t_startROM = tic;
            energyError = zeros(1,size(omega_ROM,2));
            L2Error = zeros(1,size(omega_ROM,2));
            H1Error = zeros(1,size(omega_ROM,2));
            H1sError = zeros(1,size(omega_ROM,2));
            surfaceError = zeros(1,size(omega_ROM,2));
            for i_f = 1:numel(omega_ROM)
%             parfor i_f = 1:numel(omega_ROM)
                varCol_temp = varCol;
                varCol_temp{1}.omega = omega_ROM(i_f);
                k = omega_ROM(i_f)/varCol_temp{1}.c_f;
                varCol_temp{1}.k = k;
                varCol_temp{1}.f = omega_ROM(i_f)/(2*pi);
                Uc = cell(1);
                Uc{1} = U_fluid_oArr(:,i_f);
                varCol_temp = getAnalyticSolutions(varCol_temp);
                [L2Error(i_f), H1Error(i_f), H1sError(i_f), energyError(i_f), surfaceError(i_f)] ...
                                = calculateErrors(task, varCol_temp, Uc, runTasksInParallel, stringShift);
            end
            switch basisROM
                case {'Taylor','Pade'}
                    if task.calculateSurfaceError
                        [surfaceError,temp_omega_ROM] = insertNaN(omega_ROM,surfaceError,omega_P);
                    end
                    if task.calculateVolumeError
                        [energyError,temp_omega_ROM] = insertNaN(omega_ROM,energyError,omega_P);
                        L2Error  = insertNaN(omega_ROM,L2Error,omega_P);
                        H1Error  = insertNaN(omega_ROM,H1Error,omega_P);
                        H1sError = insertNaN(omega_ROM,H1sError,omega_P);
                    end
                otherwise
                    temp_omega_ROM = omega_ROM;
            end
            if task.calculateFarFieldPattern
                varCol{1}.omega = omega_ROM;
                varCol{1}.k = omega_ROM/varCol{1}.c_f;
                varCol{1}.lambda = 2*pi./varCol{1}.k;
                varCol{1}.f = omega_ROM/(2*pi);
                Uc = cell(1);
                Uc{1} = U_fluid_oArr;
                varCol = getAnalyticSolutions(varCol);
                task = calculateTS(varCol,task,Uc,runTasksInParallel,stringShift);
                fieldCell = {'p','abs_p','TS'};
                if varCol{1}.analyticSolutionExist
                    fieldCell = [fieldCell, 'p_ref','abs_p_ref','TS_ref','error_pAbs','error_p'];
                end
                for field = fieldCell
                    tasks(i_task,taskROM,i_b).task.results.(field{1}) = task.results.(field{1});
                end
            end
            tasks(i_task,taskROM,i_b).task.results.energyError = energyError;
            tasks(i_task,taskROM,i_b).task.results.L2Error = L2Error;
            tasks(i_task,taskROM,i_b).task.results.H1Error = H1Error;
            tasks(i_task,taskROM,i_b).task.results.H1sError = H1sError;
            tasks(i_task,taskROM,i_b).task.results.surfaceError = surfaceError;
            tasks(i_task,taskROM,i_b).task.varCol = extractVarColFields(task,varCol);
            tasks(i_task,taskROM,i_b).task.varCol{1}.omega_ROM = temp_omega_ROM;
            tasks(i_task,taskROM,i_b).task.varCol{1}.f_ROM = temp_omega_ROM/(2*pi);
            tasks(i_task,taskROM,i_b).task.noVecs = noVecs;
            tasks(i_task,taskROM,i_b).task.basisROM = basisROM;
            
            fprintf('using %12f seconds.', toc(t_startROM))
        end
    end
end

function [y,x] = insertNaN(x,y,a)
inter = (a(2:end)+a(1:end-1))/2;
for i = 1:numel(a)-1
    [x,I] = sort([x, inter(i)]);
    y = [y, NaN];
    y = y(I);
end

