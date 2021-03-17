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
misc.applyLoads = {'pointPulsation','pointCharge'};
% misc.applyLoads = {'pointCharge'};
misc.applyLoads = {'pointPulsation'};
for i = 1:numel(misc.applyLoads)
    misc.applyLoad = misc.applyLoads{i};
    % misc.applyLoad = 'pointCharge';
    % misc.applyLoad = 'planeWave';
    % BCs = {'SHBC','SSBC'};
    % BCs = {'SHBC'};

    warning('off','NURBS:weights')

    prePlot.abortAfterPlotting  = 1;       % Abort simulation after pre plotting
    prePlot.plot3Dgeometry = 0;
    prePlot.plot2Dgeometry = 1;
    prePlot.plotControlPolygon = 0;       % Plot the control polygon for the NURBS mesh
    % prePlot.colorFun = @(v) abs(norm2(v)-1);
    prePlot.resolution = [20,20,0];

    postPlot(1).yname        	= 'TS';
    postPlot(1).plotResults  	= true;
    postPlot(1).printResults 	= false;
    postPlot(1).axisType      	= 'plot';
    postPlot(1).lineStyle    	= '-';
    postPlot(1).xScale       	= 1;
    postPlot(1).legendEntries 	= {};
    postPlot(1).subFolderName 	= '';
    postPlot(1).fileDataHeaderX	= [];
    postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_();

    alpha_s = 0;
    beta_s = 0;   
    alpha = alpha_s;                            % Aspect angles of observation points
    beta = beta_s;   
    if strcmp(misc.applyLoad,'pointPulsation')
        err.calculateSurfaceError = 1;
        calculateVolumeError  = 0;
        BCs = {'NBC'};
        postPlot(2) = postPlot(1);
        postPlot(2).yname        	= 'surfaceError';
        postPlot(2).axisType      	= 'semilogy';
        postPlot(2).addCommands   	= [];
    else
        BCs = {'SSBC'};
        err.calculateSurfaceError = 0;
        calculateVolumeError  = 0;
        postPlot = postPlot(1);
    end

    for BC = BCs
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
        varCol{1}.meshFile = 'createNURBSmesh_M3';
    %     k = linspace(2.5, 20, 5)/varCol{1}.R1;
%         k = linspace(0.5, 4.29, 10);
        k = linspace(0.5, 4.29, 10)/10;

        k_ROM = k(1):0.005:k(end);
    %     k_ROM = k(1):0.05:k(end);
        c_f = varCol{1}.c_f;
        omega_ROM = k_ROM*c_f;
        f = k*c_f/(2*pi);
        omega = 2*pi*f;
        
        r_a = 1.25*varCol{1}.R1;
        t_PML = 0.25*varCol{1}.R1;
        gamma = 7;
        gamma = 7*t_PML^2;
        gamma = linspace(58,79,20);
        gamma = 70.1579;
        n = 2;
        gamma = -log(1e9*eps)*(n+1)/(k(end)*t_PML);
        sigmaType = 2;

        %% Settings for the PML (perfectly matched layers)
        pml.eps = 1e9*eps;      % choosing eps = eps yields machine precicion at Gamma_b, but requires more "radial" elements in the PML to resolve the rapid decay function
        pml.sigmaType = 2;   	% sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n
        pml.t = NaN;     	    % thickness of PML
        pml.n = 2;            	% polynomial order
        pml.dirichlet = false;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
        pml.C = -log(1e9*eps)*(n+1)/(k(end)*t_PML);       
        
        postPlot(1).xScale = varCol{1}.R1;

        Xi = [0,0,0,1,1,2,2,3,3,3]/3;
        extraGP = [7,0,0];    % extra quadrature points
        varCol{1}.refinement = @(M) [0, 2^(M-1)-1, max(2^(M-1)/8-1,0), max(round(2^(M-1)-1),0)];
        if noDomains > 1
            varCol{2}.refinement = @(M,t,t_fluid) [0, 2^(M-1)-1, max(round(t/t_fluid)*2^(M-1),0)];
        end

        sdfmsc = 3*(varCol{1}.R1+varCol{1}.L+varCol{1}.R2)/1.36;
        r_s = sdfmsc - varCol{1}.L/2;
        c_z = 45;
    %     c_xy = 15;
    %     c_z = 44.45920956623927;
        c_xy = 11.439247597213537;

        basisROMcell = {'Pade','Taylor','DGP','Hermite','Bernstein'};  % do not put basisROMcell in loopParameters (this is done automatically)
        basisROMcell = {'DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
        if strcmp(misc.method,'PML')
            formulation = {'GSB'};
        else
            formulation = {'BGC'};
        end
        noVecsArr = 32;
        degree = 2;
        M = 5; % 5
        N = 20; % 9
        useROM = true;
        p_ie = 4;
        s_ie = 2;
        IElocSup = 1;        % Toggle usage of radial shape functions in IE with local support

        storeFullVarCol = false;
        if strcmp(misc.scatteringCase, 'Sweep')
            loopParameters = {'M','misc.method','BC','misc.applyLoad'};
        else
            loopParameters = {'M','misc.method','BC','misc.applyLoad','f'};
        end
%         collectIntoTasks

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     misc.scatteringCase = 'Sweep';
        misc.scatteringCase = 'BI';
        if strcmp(misc.scatteringCase,'Sweep')
            postPlot(1).noXLoopPrms   	= 0;
            if strcmp(misc.applyLoad,'pointPulsation')
                postPlot(2).noXLoopPrms = 0;
            end
            loopParameters = {'M','misc.method','BC','misc.applyLoad'};
        else
            loopParameters = {'M','misc.method','BC','misc.applyLoad','f'};
            postPlot(1).noXLoopPrms   	= 1;
            postPlot(1).xLoopName   	= 'f';
            if strcmp(misc.applyLoad,'pointPulsation')
                postPlot(2).noXLoopPrms = 1;
                postPlot(2).xLoopName = 'f';
            end
        end
        postPlot(1).xname = postPlot(1).xname(1);
        if strcmp(misc.applyLoad,'pointPulsation')
            postPlot(2).xname = postPlot(2).xname(1);
        end
        omega = omega_ROM;
        useROM = false;
        if 0 %strcmp(misc.scatteringCase, 'BI')
            para.plotResultsInParaview	 = true;	% Only if misc.scatteringCase == 'Bi'
            para.extraXiPts              = '0';  % Extra visualization points in the xi-direction per element
            para.extraEtaPts             = 'round(20/2^(M-1))';  % Extra visualization points in the eta-direction per element
            para.extraZetaPts            = 'round(1/2^(M-1))';   % Extra visualization points in the zeta-direction per element
            omega = omega_ROM(end);
        else
            omega = omega_ROM;
            omega = omega_ROM(end);
        end
        f = omega/(2*pi);
        if strcmp(misc.method,'PML')
            formulation = {'GSB'};
        else
            formulation = {'BGU'};
        end
        IElocSup = 0;        % Toggle usage of radial shape functions in IE with local support
        N = 5;
        collectIntoTasks
        
        if strcmp(misc.method,'PML')
            formulation = {'GSB'};
        else
            formulation = {'BGC'};
        end
        N = 50;
        p_ie = 4;
        s_ie = 2;
        IElocSup = 1;        % Toggle usage of radial shape functions in IE with local support
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

