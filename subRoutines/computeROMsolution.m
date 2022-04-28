function tasks = computeROMsolution(tasks,i_task,basisROMcell,omega_ROM,noVecsArr,printLog)
% basisROM = 'Taylor';
% basisROM = 'Pade';
% basisROM = 'Lagrange';
% basisROM = 'Splines';
% basisROM = 'Fourier';
% basisROM = 'Bernstein';

U_sweep = tasks(i_task).task.U_sweep;
noDofs = size(U_sweep{1},1);
A0 = tasks(i_task).task.varCol{1}.A0;
A1 = tasks(i_task).task.varCol{1}.A1;
A2 = tasks(i_task).task.varCol{1}.A2;
tasks(i_task).task.varCol = rmfields(tasks(i_task).task.varCol,{'A0','A1','A2','U'});
tasks(i_task).task = rmfield(tasks(i_task).task,'U_sweep');
task = tasks(i_task).task;
omega_P = task.misc.omega;
P = numel(omega_P);
omega_start = omega_P(1);
omega_end = omega_P(end);
stringShift = 40;
tasks(i_task).task = rmfield(tasks(i_task).task,'varCol');
noDomains = numel(task.varCol);

% noVec = size(U_P{1},2);
for i_b = 1:numel(basisROMcell)
    basisROM = basisROMcell{i_b};
    for taskROM = 1:numel(noVecsArr)
        noVecs = noVecsArr(taskROM);
        if printLog
            fprintf(['\n%-' num2str(stringShift) 's'], 'Computing basis for ROM ... ')
        end
        t_startROM = tic;
        switch basisROM
            case 'DGP'
                V = zeros(noDofs,P*noVecs);
                counter = 1;
                for i = 1:noVecs
                    for j = 1:P
                        V(:,counter) = U_sweep{j}(:,i);
                        counter = counter + 1;
                    end
                end
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
                        f = @(x,n) mp(U_sweep{i}(:,n+1));
                        [p{i},q{i}] = pade(f,mp(omega_P(i)),n,m); 
                    else
                        f = @(x,n) U_sweep{i}(:,n+1);
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
                        b(counter,:) = U_sweep{j}(:,i).';
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
                        b(counter,:) = U_sweep{j}(:,i).';
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
                        b(counter,:) = U_sweep{j}(:,i).';
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
                        b(counter,:) = U_sweep{j}(:,i).';
                        counter = counter + 1;
                    end
                end
                cond(V)
                a = V\b;
        end  
        if printLog
            fprintf('using %12f seconds.', toc(t_startROM))
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
        omega_ROM = double(omega_ROM);
        energyError = zeros(1,size(omega_ROM,2));
        L2Error = zeros(1,size(omega_ROM,2));
        H1Error = zeros(1,size(omega_ROM,2));
        H1sError = zeros(1,size(omega_ROM,2));
        surfaceError = zeros(1,size(omega_ROM,2));
        p_h = zeros(1,size(omega_ROM,2));
        for i_f = 1:numel(omega_ROM)
            omega = omega_ROM(i_f);
            task.misc.omega = omega;
            if printLog
                fprintf(['\n%-' num2str(stringShift) 's'], ['Computing ROM solution (' num2str(i_f) '/' num2str(numel(omega_ROM)) ')... '])
            end
            t_startROM = tic;
            task = getAnalyticSolutions(task);
            switch basisROM
                case 'DGP'
                    task.noRHSs = 1;
                    for i_domain = 1:min(noDomains,2)
                        task.rom.useROM = false;
                        task = applyNeumannCondition(task,i_domain);
                        task.rom.useROM = true;
                    end
                    [task,FF] = collectMatrices(task);
                    FFm = V'*FF;  

                    Am = A0_am + omega*A1_am + omega^2*A2_am;
                    Pinv = spdiags(1./diag(Am),0,size(Am,1),size(Am,2));
                    task.UU = V*(Pinv*((Am*Pinv)\FFm));
                case 'Hermite'
                    task.UU = zeros(noDofs,1);
                    counter = 1;
                    for i = 1:P
                        for n = 1:noVecs
                            task.UU = task.UU + U_sweep{i}(:,n)*Y(counter,:);
                            counter = counter + 1;
                        end
                    end
                case 'Pade'
                    if useHP
                        task.UU = double(interPade(mp(omega),mp(omega_P),p,q));
                    else
                        task.UU = interPade(omega,omega_P,p,q);
                    end
                case 'Bernstein'
                    B = bernsteinBasis(double((omega-omega_start)/(omega_end - omega_start)),double(p_ROM),0);
                    task.UU = (B*a).';
                case 'Taylor'
                    task.UU = interTaylor(omega,omega_P,U_sweep,noVecs-1);
            end
            task = computeDerivedQuantities(task);
            task = postProcessSolution(task);
            if printLog
                fprintf('using %12f seconds.', toc(t_startROM))
            end

            if task.err.calculateSurfaceError || task.err.calculateVolumeError
                [L2Error(i_f), H1Error(i_f), H1sError(i_f), energyError(i_f), surfaceError(i_f)] = calculateErrors(task, 1, stringShift);
            end

            if task.ffp.calculateFarFieldPattern
                task = calculateTS(task,printLog,stringShift);
                p_h(i_f) = task.ffp.p_h;
            end
        end
        if task.ffp.calculateFarFieldPattern
            task.misc.omega = omega_ROM;
            task = getAnalyticSolutions(task);
            task = computeDerivedFFPquantities(task,p_h);
            fieldCell = {'p','p_Re','p_Im','abs_p','TS'};
            if task.analyticSolutionExist
                fieldCell = [fieldCell, 'p_ref', 'p_Re_ref', 'p_Im_ref','abs_p_ref','TS_ref','error_pAbs','error_p'];
            end
            for field = fieldCell
                switch basisROM
                    case {'Taylor','Pade'}
                        [tasks(i_task,taskROM,i_b).task.results.(field{1}),temp_omega_ROM] = insertNaN(omega_ROM,task.results.(field{1}),omega_P);
                    otherwise
                        tasks(i_task,taskROM,i_b).task.results.(field{1}) = task.results.(field{1});
                        temp_omega_ROM = omega_ROM;
                end
            end
        end
        if task.err.calculateSurfaceError || task.err.calculateVolumeError
            switch basisROM
                case {'Taylor','Pade'}
                    if task.err.calculateSurfaceError
                        [surfaceError,temp_omega_ROM] = insertNaN(omega_ROM,surfaceError,omega_P);
                    end
                    if task.err.calculateVolumeError
                        [energyError,temp_omega_ROM] = insertNaN(omega_ROM,energyError,omega_P);
                        L2Error  = insertNaN(omega_ROM,L2Error,omega_P);
                        H1Error  = insertNaN(omega_ROM,H1Error,omega_P);
                        H1sError = insertNaN(omega_ROM,H1sError,omega_P);
                    end
                otherwise
                    temp_omega_ROM = omega_ROM;
            end
            tasks(i_task,taskROM,i_b).task.results.energyError = energyError;
            tasks(i_task,taskROM,i_b).task.results.L2Error = L2Error;
            tasks(i_task,taskROM,i_b).task.results.H1Error = H1Error;
            tasks(i_task,taskROM,i_b).task.results.H1sError = H1sError;
            tasks(i_task,taskROM,i_b).task.results.surfaceError = surfaceError;
        end

        tasks(i_task,taskROM,i_b).task = rmfield(task,'varCol');
        tasks(i_task,taskROM,i_b).task.varCol{1}.omega_ROM = temp_omega_ROM;
        tasks(i_task,taskROM,i_b).task.varCol{1}.f_ROM = temp_omega_ROM/(2*pi);
        tasks(i_task,taskROM,i_b).task.varCol{1}.k_ROM = temp_omega_ROM/task.varCol{1}.c_f;
        tasks(i_task,taskROM,i_b).task.rom.noVecs = noVecs;
        tasks(i_task,taskROM,i_b).task.rom.basisROM = basisROM;
        fieldCell = {'misc','msh','sol','err','ffp','iem','pml','bem','rt','mfs','dofs','noRHSs','analyticSolutionExist','isSphericalShell','totNoElems','ffp'};
        for field = fieldCell
            tasks(i_task,taskROM,i_b).task.(field{1}) = task.(field{1});
        end
            
        tasks(i_task,taskROM,i_b).task.msh = task.msh;
        tasks(i_task,taskROM,i_b).task.misc = task.misc;
        tasks(i_task,taskROM,i_b).task.sol = task.sol;
        tasks(i_task,taskROM,i_b).task.err = task.err;
        tasks(i_task,taskROM,i_b).task.ffp = task.ffp;
        tasks(i_task,taskROM,i_b).task.iem = task.iem;
        tasks(i_task,taskROM,i_b).task.pml = task.pml;
        tasks(i_task,taskROM,i_b).task.bem = task.bem;
        tasks(i_task,taskROM,i_b).task.rt = task.rt;
        tasks(i_task,taskROM,i_b).task.mfs = task.mfs;
        tasks(i_task,taskROM,i_b).task.dofs = task.dofs;
    end
end

function [y,x] = insertNaN(x,y,a)
inter = (a(2:end)+a(1:end-1))/2;
for i = 1:numel(a)-1
    [x,I] = sort([x, inter(i)]);
    y = [y, NaN];
    y = y(I);
end

