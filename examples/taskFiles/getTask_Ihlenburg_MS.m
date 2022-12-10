function studies = getTask_Ihlenburg_MS()
% This study is based on Ihlenburg1998fea (Figure 7)
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271
% Experimental data at scale 1:50
counter = 1;
studies = cell(0,1);
getDefaultTaskValues

saveStudies = false;
noCoresToUse = 32;         % Number of processors for parallel computations (Inf uses all available cores)
% noCoresToUse = 4;         % Number of processors for parallel computations (Inf uses all available cores)


%% IE simulation
misc.model = 'IMS';  % Spherical shell

misc.coreMethod = {'IGA','hp_FEM'};
misc.coreMethod = {'IGA'};
% misc.method = {'IE'};
misc.method = {'PML'};
applyLoads = {'pointPulsation','pointCharge'};
applyLoads = {'pointCharge'};
% applyLoads = {'pointPulsation'};

msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
msh.refineThetaOnly = 1;

misc.checkNURBSweightsCompatibility = false;
misc.preProcessOnly = 0;
warning('off','NURBS:weights')

prePlot.plot3Dgeometry = 0;
prePlot.abortAfterPlotting = 1;       % Abort simulation after pre plotting
prePlot.plotControlPolygon = 0;       % Plot the control polygon for the NURBS mesh
% prePlot.colorFun = @(v) abs(norm2(v)-1);
%     prePlot.resolution = [20,20,0];
prePlot.resolution = [400,100,0];
prePlot.resolution = [120,40,0];
prePlot.plotFullDomain      = 0;
prePlot.plotSubsets         = {'xy'}; 
% prePlot.plotSubsets         = {'innerCoupling','outerCoupling','inner','outer'}; 
prePlot.view                = [0,90];     % Set view angle [azimuth,elevation]
prePlot.useCamlight         = false;
% prePlot.format              = '-pdf';      % Use '-png' or '-pdf' (for vector graphics)
prePlot.pngResolution       = '-r600';      % Use '-png' or '-pdf' (for vector graphics)

postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType      	= 'plot';
postPlot(1).lineStyle    	= '-';
postPlot(1).legendEntries 	= {};
postPlot(1).subFolderName 	= '';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(study);

ffp.alpha_s = 0;
ffp.beta_s = 0;   
ffp.alpha = ffp.alpha_s;                            % Aspect angles of observation points
ffp.beta = ffp.beta_s;  
msh.meshFile = 'createNURBSmesh_M3'; 
if msh.refineThetaOnly
    ffp.extraGP = [100,0,0];                      % Extra Gauss points used for the integration routine
end
for i = 1:numel(applyLoads)
    misc.applyLoad = applyLoads{i};
    % BCs = {'SHBC','SSBC'};
    % BCs = {'SHBC'};
    if strcmp(misc.applyLoad,'pointPulsation')
        err.calculateSurfaceError = 1;
        err.calculateVolumeError  = 0;
        BCs = {'NBC'};
        postPlot(2) = postPlot(1);
        postPlot(2).yname        	= 'surfaceError';
        postPlot(2).axisType      	= 'semilogy';
        postPlot(2).addCommands   	= [];
    else
        BCs = {'SSBC'};
%         BCs = {'SHBC'};
        err.calculateSurfaceError = 0;
        err.calculateVolumeError  = 0;
        postPlot = postPlot(1);
    end

    for BC = BCs
        misc.BC = BC;
        misc.scatteringCase = 'Sweep';
        postPlot(1).noXLoopPrms = 0;
        postPlot(1).xname = 'varCol{1}.k';
        if strcmp(misc.applyLoad,'pointPulsation')
            postPlot(2).xname = 'varCol{1}.k';
        end
        switch BC{1}
            case {'SHBC','NBC'}
                noDomains = 1;
            case 'SSBC'
                noDomains = 2;
        end
        varCol{1} = struct('media', 'fluid', ...
                           't', [0.1143,0.04064], ...
                           'R1', 4.6863, ...
                           'R2', 4.6863, ...
                           'L', 67.9958, ...
                           'c_f', 1482, ...
                           'rho', 1000);
        eta = 0.001; % hysteretic loss factor for hysteresis damping
        E = 2.0e11;
%         E = E*(1-1i*eta);
        varCol{2} = struct('media', 'solid', ...
                           'E', E, ...
                           'nu', 0.29, ...
                           'rho', 7908.5);
    %     varCol{1} = struct('media', 'fluid', ...
    %                        't', [0.121979301545535,0.056298757842532], ...
    %                        'R1', 4.689073447525605, ...
    %                        'R2', 4.689073447525605, ...
    %                        'L', 33.96095231541791*2, ...
    %                        'c_f', 1482, ...
    %                        'rho', 1000);
    %     varCol{2} = struct('media', 'solid', ...
    %                        'E', 2.0e11, ...
    %                        'nu', 0.29, ...
    %                        'rho', 7908.5);

        varCol = varCol(1:noDomains);
        n = 20;
%         n = 3;
        eqDistr = 1:n;
%         k = linspace(2.5, 20, n)/varCol{1}.R1;
        k_P = linspace(0.5, 4.29, n);

        k = k_P(1):0.005:k_P(end);
%         k = k(1):0.05:k(end);
%         k = k(1):0.5:k(end);
        k = sort(unique([k,k_P]));
%         k = k;
%         k = k(1)+(k(end)-k(1))*(1-cos((2*eqDistr-1)/2/n*pi))/2; % Chebyshev nodes
        c_f = varCol{1}.c_f;
        misc.omega = k*c_f;
        f_P = k_P*c_f/(2*pi);
        rom.omega = 2*pi*f_P;
        msh.explodeNURBS = 0;   % Create patches from all C^0 interfaces
        
        %% Settings for the PML (perfectly matched layers)
        refLength = varCol{1}.R1*pi/2;
        pml.sigmaType = 3;   	% sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n
        pml.t = 0.25*varCol{1}.R1; % thickness of PML
        pml.n = 1;            	% polynomial order
        pml.dirichlet = true;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
        pml.refinement = @(M) max(round((2^(M-1)-1)*pml.t/refLength),2);   
        pml.gamma = 1/(k_P(1)*pml.t);     
        
        varCol{1}.refinement = @(M) [0, 2^(M-1)-1, max(2^(M-4)-1,2)];
        if noDomains > 1
            varCol{2}.refinement = @(M,t,t_fluid) [0, 2^(M-1)-1, max(ceil(t/t_fluid)*2^(M-1),0)];
        end
        scale = (varCol{1}.R1+varCol{1}.L+varCol{1}.R2)/1.36; % Scaling between experimental and theoretical setup
        msh.x_0 = [-varCol{1}.L/2,0,0];       % Move the origin to the center of the model
        misc.r_s = 3*scale;                   % Distance from the center of the model to the point charge
        msh.c_z = 45;
    %     c_xy = 15;
    %     c_z = 44.45920956623927;
        msh.c_xy = 11.439247597213537;

%         rom.basisROM = {'Pade','Taylor','DGP','Hermite','Bernstein'};  % do not put basisROMcell in loopParameters (this is done automatically)
        rom.basisROM = {'DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
        rom.adaptiveROM  = true;
        rom.computeROMresidualFine = true;
        rom.computeROMerror = true;
        if strcmp(misc.method,'PML')
            misc.formulation = {'GSB'};
        else
            misc.formulation = {'BGC'};
        end
        msh.degree = 3:4;
%         msh.degree = 4; % 4
        msh.M = 7:8; % 7
%         msh.M = 4; % 7
        misc.symmetric = 0;
        
        misc.extraGP = [9-msh.degree(1),0,0];    % extra quadrature points
        
        rom.useROM = true;
        rom.adaptiveROM = true;
%         rom.noVecs = [4,8,16,32];
        rom.noVecs = 16;

        misc.r_a = 1.25*varCol{1}.R1;
        postPlot(1).xScale = varCol{1}.R1;
        postPlot(1).xLabel = '$$ka$$';
        if strcmp(misc.applyLoad,'pointPulsation')
            postPlot(2).xScale = varCol{1}.R1;
            postPlot(2).xLabel = '$$ka$$';
        end
        
        misc.storeFullVarCol = false;
        if strcmp(misc.scatteringCase, 'Sweep')
            loopParameters = {'msh.M','msh.degree','misc.method','misc.BC','misc.applyLoad','rom.noVecs','rom.basisROM'};
        else
            loopParameters = {'msh.M','msh.degree','misc.method','misc.BC','misc.applyLoad','rom.noVecs','rom.basisROM','misc.omega'};
        end
        collectIntoTasks

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        misc.scatteringCase = 'BI';
        if strcmp(misc.scatteringCase,'Sweep')
            postPlot(1).noXLoopPrms   	= 0;
            if strcmp(misc.applyLoad,'pointPulsation')
                postPlot(2).noXLoopPrms = 0;
            end
            loopParameters = {'msh.M','msh.degree','misc.method','misc.BC','misc.applyLoad'};
        else
            loopParameters = {'msh.M','msh.degree','misc.method','misc.BC','misc.applyLoad','misc.omega'};
            postPlot(1).noXLoopPrms   	= 1;
            postPlot(1).xLoopName   	= 'misc.omega';
            if strcmp(misc.applyLoad,'pointPulsation')
                postPlot(2).noXLoopPrms = 1;
                postPlot(2).xLoopName = 'misc.omega';
            end
        end
        postPlot(1).xname = 'varCol{1}.k';
        if strcmp(misc.applyLoad,'pointPulsation')
            postPlot(2).xname = 'varCol{1}.k';
        end
        rom.useROM = false;
        if 1 %strcmp(misc.scatteringCase, 'BI')
            para.plotResultsInParaview	 = true;	% Only if misc.scatteringCase == 'Bi'
            para.extraXiPts              = '60';  % Extra visualization points in the xi-direction per element
            para.extraEtaPts             = '1';  % Extra visualization points in the eta-direction per element
            para.extraZetaPts            = '1';   % Extra visualization points in the zeta-direction per element
            para.plotTimeOscillation     = 1;
            misc.omega = misc.omega(end);
        else
%             misc.omega = misc.omega(end);
        end
        if strcmp(misc.method,'PML')
            misc.formulation = {'GSB'};
        else
            misc.formulation = {'BGU'};
        end
        iem.IElocSup = 0;        % Toggle usage of radial shape functions in IE with local support
        iem.N = 5;
%         collectIntoTasks
    end
end

function addCommands_(study)
T = readtable('miscellaneous/refSolutions/IMS.csv','FileType','text', 'HeaderLines',0);
x = T.Var1;
y = T.Var2;
plot(x,y,'DisplayName','Experiment')
legend('off');
legend('show');
hold on

for i_task = 1:numel(study.tasks)
    task = study.tasks(i_task).task;
    c_f = task.varCol{1}.c_f;

    if isfield(task.rom,'history')
        history = task.rom.history;
        a = task.varCol{1}.R1;
        options.xlabel = 'ka';
        options.ylabel = 'residual';
        for i = 1:numel(history)
            figure(20+i)
            omega = task.rom.history(i).omega;
            residual = task.rom.history(i).residual;
            hold on
        
            omega_P = history(i).omega_P;
            J_P = history(i).J_P;
            semilogy(a*omega_P/c_f,1e-15*ones(size(omega_P)),'o','color','magenta','DisplayName','Interpolation points')
            options.x = a*omega_P.'/c_f;
            options.y = 1e-15*ones(size(omega_P)).';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_InterpolationPoints_iter' num2str(i)], options)
            for j = 1:numel(omega_P)
                text(a*(omega_P(j)+omega_P(end)/200)/c_f,1e-15,num2str(J_P(j)),'color','magenta')
            end
            options.x = a*omega_P.'/c_f;
            options.y = J_P.';
            options.ylabel = 'J_P';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_noDerivatives_iter' num2str(i)], options)
    
            semilogy(a*[omega_P(1),omega_P(end)]/c_f,task.rom.tolerance*ones(1,2),'red','DisplayName','Tolerance')
            options.x = a*[omega_P(1),omega_P(end)].'/c_f;
            options.y = task.rom.tolerance*ones(2,1);
            options.ylabel = 'residual';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_tolerance'], options)
    
            semilogy(a*omega/c_f,residual,'*','color','black','DisplayName','Residual check points')
            options.x = a*omega.'/c_f;
            options.y = residual.';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_residual_iter' num2str(i)], options)
    
            omega_T_new = history(i).omega_T_new;
            [~,idx] = ismember(omega,omega_T_new);
            semilogy(a*omega_T_new/c_f,residual(logical(idx)),'o','color','cyan','DisplayName','New interpolation point')
    
            options.x = a*omega_T_new.'/c_f;
            options.y = residual(logical(idx)).';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_NewInterpolationPoints_iter' num2str(i)], options)
    
            ylim([5e-16,10])
            set(gca,'yscale','log')
            legend show
            savefig([task.resultsFolder, '/', task.saveName '_adaptiveROM_iter' num2str(i)])
        end
        legend off
        if isfield(task.rom.history(end),'residualFine')
            residual = task.rom.history(end).residualFine;
            omegaFine = task.rom.history(end).omegaFine;
            semilogy(a*omegaFine/c_f,residual,'blue','DisplayName','Final relative residual')
        
            options.x = a*omegaFine.'/c_f;
            options.y = residual.';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_FinalRelativeResidual'], options)
        end
        if isfield(task.rom.history(end),'relError')
            relError = task.rom.history(end).relError;
            semilogy(a*omegaFine/c_f,relError,'green','DisplayName','Final relative error')
        
            options.x = a*omegaFine.'/c_f;
            options.y = relError.';
            options.xlabel = 'k';
            options.ylabel = 'error';
            printResultsToFile([task.resultsFolder, '/', task.saveName '_FinalRelativeError'], options)
        end
        set(gca,'yscale','log')
        ylim([5e-16,10])
        ylabel('Relative error/residual')
        xlabel('Wavenumber')
        legend show
        savefig([task.resultsFolder, '/', task.saveName '_adaptiveROM_iter' num2str(i)])
    end
end


