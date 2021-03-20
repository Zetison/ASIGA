function studies = getTask_Ihlenburg_MS()
% This study is based on Ihlenburg1998fea (Figure 7)
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271
% Experimental data at scale 1:50
counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
misc.model = 'IMS';  % Spherical shell

misc.coreMethod = {'IGA','hp_FEM'};
misc.coreMethod = {'IGA'};
% misc.method = {'IE'};
misc.method = {'PML'};
applyLoads = {'pointPulsation','pointCharge'};
applyLoads = {'pointCharge'};
% applyLoads = {'pointPulsation'};
for i = 1:numel(applyLoads)
    misc.applyLoad = applyLoads{i};
    % BCs = {'SHBC','SSBC'};
    % BCs = {'SHBC'};

    warning('off','NURBS:weights')

    prePlot.abortAfterPlotting  = 1;       % Abort simulation after pre plotting
    prePlot.plot3Dgeometry = 0;
    prePlot.plot2Dgeometry = 0;
    prePlot.plotControlPolygon = 0;       % Plot the control polygon for the NURBS mesh
    % prePlot.colorFun = @(v) abs(norm2(v)-1);
    prePlot.resolution = [20,20,0];
%     prePlot.resolution = [0,0,0];

    postPlot(1).yname        	= 'TS';
    postPlot(1).plotResults  	= true;
    postPlot(1).printResults 	= true;
    postPlot(1).axisType      	= 'plot';
    postPlot(1).lineStyle    	= '-';
    postPlot(1).xScale       	= 1;
    postPlot(1).legendEntries 	= {};
    postPlot(1).subFolderName 	= '';
    postPlot(1).fileDataHeaderX	= [];
    postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_();

    ffp.alpha_s = 0;
    ffp.beta_s = 0;   
    ffp.alpha = ffp.alpha_s;                            % Aspect angles of observation points
    ffp.beta = ffp.beta_s;   
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
        postPlot(1).xname = 'k_ROM';
        if strcmp(misc.applyLoad,'pointPulsation')
            postPlot(2).xname = 'k_ROM';
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
        varCol{2} = struct('media', 'solid', ...
                           'E', 2.0e11, ...
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
        msh.meshFile = 'createNURBSmesh_M3';
    %     k = linspace(2.5, 20, 5)/varCol{1}.R1;
        k = linspace(0.5, 4.29, 20);
%         k = linspace(0.5, 4.29, 3);
%         k = linspace(0.5, 4.29, 10)/10;

        k_ROM = k(1):0.005:k(end);
%         k_ROM = k(1):0.05:k(end);
%         k_ROM = k(1):0.5:k(end);
        k_ROM = sort(unique([k_ROM,k]));
%         k_ROM = k;
        c_f = varCol{1}.c_f;
        omega_ROM = k_ROM*c_f;
        f = k*c_f/(2*pi);
        misc.omega = 2*pi*f;
        
        %% Settings for the PML (perfectly matched layers)
        pml.eps = 1e0*eps;      % choosing eps = eps yields machine precicion at Gamma_b, but requires more "radial" elements in the PML to resolve the rapid decay function
        pml.sigmaType = 2;   	% sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n
        pml.t = 0.25*varCol{1}.R1; % thickness of PML
        pml.n = 2;            	% polynomial order
        pml.dirichlet = true;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
        pml.C = -log(1e9*eps)*(pml.n+1)/(k(1)*pml.t);     
        
        postPlot(1).xScale = varCol{1}.R1;

        msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
        msh.refineThetaOnly = true;
        varCol{1}.refinement = @(M) [0, 2^(M-1)-1, max(2^(M-1)/8-1,0), max(round(2^(M-1)-1),0)];
        if noDomains > 1
            varCol{2}.refinement = @(M,t,t_fluid) [0, 2^(M-1)-1, max(round(t/t_fluid)*2^(M-1),0)];
        end
        scale = (varCol{1}.R1+varCol{1}.L+varCol{1}.R2)/1.36; % Scaling between experimental and theoretical setup
        msh.x_0 = [-varCol{1}.L/2,0,0];       % Move the origin to the center of the model
        misc.r_s = 3*scale;                   % Distance from the center of the model to the point charge
        msh.c_z = 45;
    %     c_xy = 15;
    %     c_z = 44.45920956623927;
        msh.c_xy = 11.439247597213537;

        basisROMcell = {'Pade','Taylor','DGP','Hermite','Bernstein'};  % do not put basisROMcell in loopParameters (this is done automatically)
        basisROMcell = {'DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
        if strcmp(misc.method,'PML')
            misc.formulation = {'GSB'};
        else
            misc.formulation = {'BGC'};
        end
        msh.degree = 4;
        msh.M = 5; % 5
        
        misc.extraGP = [9-msh.degree,0,0];    % extra quadrature points
        
        iem.N = 20; % 9
        iem.p_ie = 4;
        iem.s_ie = 2;
        iem.IElocSup = 1;        % Toggle usage of radial shape functions in IE with local support
        
        rom.useROM = true;
        rom.noVecsArr = 32;
%         rom.noVecsArr = 1;

        misc.r_a = 6.3158;
        
        misc.storeFullVarCol = false;
        if strcmp(misc.scatteringCase, 'Sweep')
            loopParameters = {'msh.M','misc.method','misc.BC','misc.applyLoad'};
        else
            loopParameters = {'msh.M','misc.method','misc.BC','misc.applyLoad','misc.omega'};
        end
%         collectIntoTasks

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        misc.scatteringCase = 'BI';
        if strcmp(misc.scatteringCase,'Sweep')
            postPlot(1).noXLoopPrms   	= 0;
            if strcmp(misc.applyLoad,'pointPulsation')
                postPlot(2).noXLoopPrms = 0;
            end
            loopParameters = {'msh.M','misc.method','misc.BC','misc.applyLoad'};
        else
            loopParameters = {'msh.M','misc.method','misc.BC','misc.applyLoad','misc.omega'};
            postPlot(1).noXLoopPrms   	= 1;
            postPlot(1).xLoopName   	= 'misc.omega';
            if strcmp(misc.applyLoad,'pointPulsation')
                postPlot(2).noXLoopPrms = 1;
                postPlot(2).xLoopName = 'misc.omega';
            end
        end
        postPlot(1).xname = postPlot(1).xname(1);
        if strcmp(misc.applyLoad,'pointPulsation')
            postPlot(2).xname = postPlot(2).xname(1);
        end
        misc.omega = omega_ROM;
        rom.useROM = false;
        if 1 %strcmp(misc.scatteringCase, 'BI')
            para.plotResultsInParaview	 = true;	% Only if misc.scatteringCase == 'Bi'
            para.extraXiPts              = '60';  % Extra visualization points in the xi-direction per element
            para.extraEtaPts             = '1';  % Extra visualization points in the eta-direction per element
            para.extraZetaPts            = '1';   % Extra visualization points in the zeta-direction per element
            para.plotTimeOscillation     = 1;
            misc.omega = omega_ROM(end);
        else
            misc.omega = omega_ROM;
%             misc.omega = omega_ROM(end);
        end
        if strcmp(misc.method,'PML')
            misc.formulation = {'GSB'};
        else
            misc.formulation = {'BGU'};
        end
        iem.IElocSup = 0;        % Toggle usage of radial shape functions in IE with local support
        iem.N = 5;
%         collectIntoTasks
        
        iem.N = 50;
        iem.p_ie = 4;
        iem.s_ie = 2;
        iem.IElocSup = 1;        % Toggle usage of radial shape functions in IE with local support
    %     collectIntoTasks
    end
end
function addCommands_()
T = readtable('miscellaneous/refSolutions/IMS.csv','FileType','text', 'HeaderLines',1);
x = T.Var1;
y = T.Var2;
plot(x,y,'DisplayName','Experiment')
legend('off');
legend('show');
hold on

