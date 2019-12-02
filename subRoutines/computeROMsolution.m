
% basisROM = 'Taylor';
% basisROM = 'Pade';
% basisROM = 'Lagrange';
% basisROM = 'Splines';
% basisROM = 'Fourier';
% basisROM = 'Bernstein';
U_P = task.varCol{1}.U_sweep;
stringShift = 40;
task.varCol{1} = rmfield(task.varCol{1},'U_sweep');
task.varCol{1} = rmfield(task.varCol{1},'U');
runTasksInParallelOld = runTasksInParallel;
basisROMcell = studies(study_i).basisROMcell;
k_ROM = studies(study_i).k_ROM;
noVecsArr = studies(study_i).noVecsArr;
noDofs = size(U_P{1},1);
varCol = task.varCol;
e3Dss_options = varCol{1}.e3Dss_options;
% noVec = size(U_P{1},2);
for i_b = 1:numel(basisROMcell)
    basisROM = basisROMcell{i_b};
    k_P = task.varCol{1}.k;

    P = numel(k_P);
    k_start = k_P(1);
    k_end = k_P(end);
    for task_ROM = 1:numel(noVecsArr)
        noVecs = noVecsArr(task_ROM);

        fprintf(['\n%-' num2str(stringShift) 's'], 'Computing basis for ROM ... ')
        t_startROM = tic;
        switch basisROM
            case 'Hermite'
                mp.Digits(400);
                Y = getInterpolatingHermite(mp(k_P.'),mp(k_ROM),noVecs);
                Y = double(Y);
            case 'Pade'
                p = cell(P,1);
                q = cell(P,1);
                useHP = 0;
                if useHP
%                     mp.Digits(400);
                    mp.Digits(1000);
                end
                for i = 1:P
                    n = ceil(noVecs/2)-1;
                    m = floor(noVecs/2);
                    if useHP
                        f = @(x,n) mp(U_P{i}(:,n+1));
                        [p{i},q{i}] = pade(f,mp(k_P(i)),n,m); 
                    else
                        f = @(x,n) U_P{i}(:,n+1);
                        [p{i},q{i}] = pade(f,k_P(i),n,m); 
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
                    xi = (k_P(i)-k_start)/(k_end - k_start);
                    i1 = findKnotSpan(Vsize, p_ROM, xi, Xi);
                    ders = Bspline_basisDers3(i1, xi, p_ROM, Xi, noVecs-1);
                    for j = 2:noVecs
                        ders(j,:) = ders(j,:)/(k_end - k_start).^(j-1);
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
                k_P = convert(k_P,class(p_ROM));
                k_end = k_P(end);
                k_start = k_P(1);
                tic
                B = bernsteinBasis((k_P-k_start)/(k_end - k_start),p_ROM,noVecs-1,1/(k_end - k_start));
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
                a = double(V\mp(b));
                toc
        %         a = double(invV)*b;

            case 'Lagrange'
                Vsize = P*noVecs;
                P = noTasks;
                n = Vsize;
                temp = cos(pi/(2*n));
                ak = ((k_start+k_end)*temp-(k_end-k_start))/(2*temp);
                bk = ((k_start+k_end)*temp+(k_end-k_start))/(2*temp);
                j = 1:n;
                k_arrLagr = 1/2*(ak+bk)+1/2*(bk-ak)*cos((2*n-2*j+1)*pi/(2*n));

                dp = zeros(P,noVecs,n);
                for i = 1:n
                    dp(:,1,i) = lagrangePolynomials(k_P.',i,n,k_arrLagr);
                    dp(:,2:noVecs,i) = lagrangePolynomialsNthDeriv(k_P.',i,n,k_arrLagr,noVecs-1);
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
                            V(counter,2:end) = (-1)^((i-1)/2)*(n_arr*pi/(k_end-k_start)).^(i-1).*cos(n_arr*pi*k_P(j)/(k_end-k_start));
                        else
                            V(counter,2:end) = (-1)^(i/2)*(n_arr*pi/(k_end-k_start)).^(i-1).*sin(n_arr*pi*k_P(j)/(k_end-k_start));
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
%             k_ROM = double(unique(sort([k_P,k_ROM])));
            k_ROM = double(k_ROM);
            fprintf(['\n%-' num2str(stringShift) 's'], 'Computing ROM solution ... ')
            t_startROM = tic;
            switch basisROM
                case 'Hermite'
                    U_fluid_oArr = zeros(noDofs,numel(k_ROM));
                    counter = 1;
                    for i = 1:P
                        for n = 1:noVecs
                            U_fluid_oArr = U_fluid_oArr + U_P{i}(:,n)*Y(counter,:);
                            counter = counter + 1;
                        end
                    end
                case 'Pade'
                    if useHP
                        U_fluid_oArr = double(interPade(mp(k_ROM),mp(k_P),p,q));
                    else
                        U_fluid_oArr = interPade(k_ROM,k_P,p,q);
                    end
                case 'Bernstein'
                    B = bernsteinBasis(double((k_ROM-k_start)/(k_end - k_start)),double(p_ROM),0);
                    U_fluid_oArr = (B*a).';
                case 'Taylor'
                    U_fluid_oArr = interTaylor(k_ROM,k_P,U_P,noVecs-1);
            end
            fprintf('using %12f seconds.', toc(t_startROM))
            surfaceErrorArr = zeros(size(k_ROM));
            runTasksInParallel = 0; % to avoid printing elapsed time
            fprintf(['\n%-' num2str(stringShift) 's'], 'Computing errors for ROM sweeps ... ')
            t_startROM = tic;
            for i_k = 1:numel(k_ROM)
                k = k_ROM(i_k);
                varCol{1}.k = k;
                omega = varCol{1}.c_f*k;
                varCol{1}.f = omega/(2*pi);
                Uc{1} = U_fluid_oArr(:,i_k);
                e3Dss_options.omega = omega;
                varCol{1}.analytic = @(x)analytic_(x,e3Dss_options);
                calculateErrors
            end
            fprintf('using %12f seconds.', toc(t_startROM))
            tasks(i_task,task_ROM,i_b).task = task;
            switch basisROM
                case {'Taylor','Pade'}
                    if task.calculateSurfaceError
                        [temp_error,temp_k_ROM] = insertNaN(k_ROM,task.results.surfaceError,k_P);
                        tasks(i_task,task_ROM,i_b).task.varCol{1}.k_ROM = temp_k_ROM;
                        tasks(i_task,task_ROM,i_b).task.results.surfaceError = temp_error;
                    end
                    if task.calculateVolumeError
                        [temp_error,temp_k_ROM] = insertNaN(k_ROM,task.results.energyError,k_P);
                        tasks(i_task,task_ROM,i_b).task.varCol{1}.k_ROM = temp_k_ROM;
                        tasks(i_task,task_ROM,i_b).task.results.energyError = temp_error;
                        temp_error = insertNaN(k_ROM,task.results.L2Error,k_P);
                        tasks(i_task,task_ROM,i_b).task.results.L2Error = temp_error;
                        temp_error = insertNaN(k_ROM,task.results.H1Error,k_P);
                        tasks(i_task,task_ROM,i_b).task.results.H1Error = temp_error;
                        temp_error = insertNaN(k_ROM,task.results.H1sError,k_P);
                        tasks(i_task,task_ROM,i_b).task.results.H1sError = temp_error;
                    end
                otherwise
                    if task.calculateSurfaceError
                        tasks(i_task,task_ROM,i_b).task.results.surfaceError = task.results.surfaceError;
                    end
                    if task.calculateVolumeError
                        tasks(i_task,task_ROM,i_b).task.results.energyError = task.results.energyError;
                        tasks(i_task,task_ROM,i_b).task.results.L2Error = task.results.L2Error;
                        tasks(i_task,task_ROM,i_b).task.results.H1Error = task.results.H1Error;
                        tasks(i_task,task_ROM,i_b).task.results.H1sError = task.results.H1sError;
                        tasks(i_task,task_ROM,i_b).task.varCol{1}.k_ROM = k_ROM;
                    end
            end
            tasks(i_task,task_ROM,i_b).task.noVecs = noVecs;
            tasks(i_task,task_ROM,i_b).task.basisROM = basisROM;

%             figure(52)
%             semilogy(k_ROM,surfaceErrorArr,'DisplayName',['Bernstein-d - ' num2str(noVecs) '+' num2str(noVecs) '+' num2str(noVecs) ' vectors'])
%             hold on
%             drawnow
        end
    end
%     legend show
    % savefig('ROM')
end
runTasksInParallel = runTasksInParallelOld;

function [y,x] = insertNaN(x,y,a)
inter = (a(2:end)+a(1:end-1))/2;
for i = 1:numel(a)-1
    [x,I] = sort([x, inter(i)]);
    y = [y, NaN];
    y = y(I);
end

end

