function studies = getTask_Hetmaniuk()
% This study is based on Hetmaniuk2012raa (Figure 7)
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271

counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
hetmaniukCase = 1; % evaluating solution at boundary not implemented

% BC = 'NNBC';
misc.model = 'Hetmaniuk2012raa';  % Spherical shell
misc.coreMethod = {'IGA','C0_IGA'};
% misc.coreMethod = {'C0_IGA'};
misc.coreMethod = {'IGA'};
% misc.coreMethod = {'hp_FEM'};
% misc.applyLoad = 'pointPulsation';
% misc.applyLoad = 'pointCharge';
misc.applyLoad = 'planeWave';

% BCs = {'SHBC'};
BCs = {'SSBC'};
% BCs = {'NNBC'};
BCs = {'SHBC','SSBC'};
if strcmp(misc.applyLoad,'pointPulsation')
    BCs = {'NBC'};
end
ffp.plotFarField = ~hetmaniukCase;
% plotFarField = true;     % If false, plots the near field instead

ffp.calculateFarFieldPattern    = true;     % Calculate far field pattern
ffp.alpha_s = 0;                            % Aspect angle of incident wave
ffp.beta_s  = -pi/2;                        % Elevation angle of incident wave
ffp.alpha   = 0;                            % Aspect angles of observation points
ffp.beta = pi/2;   
ffp.r = 1;                            % radii for near-field evaluation.
ffp.splineBasedNFPcalc = true;

misc.r_s = 3;

msh.parm = 1;
err.calculateSurfaceError = 1;
err.calculateVolumeError  = 0;
misc.calculateFarFieldPattern = 1;

prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;
% prePlot.colorFun = @(v) abs(norm2(v)-1);
prePlot.resolution = [100,40,0];

misc.computeCondNumber = 0;

for i = 1:numel(BCs)
    misc.BC = BCs{i};
%     misc.method = {'IENSG'};
    misc.method = {'IE'};
%     misc.method = {'PML'};

    postPlot(1).xname           = 'k_ROM';
    postPlot(1).yname        	= 'surfaceError';
    postPlot(1).plotResults  	= true;
    postPlot(1).printResults 	= false;
    postPlot(1).axisType      	= 'semilogy';
    postPlot(1).lineStyle    	= '-';
    postPlot(1).xScale       	= 1;
    postPlot(1).yScale       	= 1;
    postPlot(1).legendEntries 	= {};
    postPlot(1).subFolderName 	= '';
    postPlot(1).fileDataHeaderX	= [];
    postPlot(1).noXLoopPrms   	= 0;
    postPlot(1).xLoopName       = NaN;

    postPlot(2) = postPlot(1);
    if hetmaniukCase
        postPlot(2).yname = 'p_Re';
    else
        postPlot(2).yname = 'TS';
    end
    postPlot(2).axisType      	= 'plot';

    postPlot(3) = postPlot(2);
    postPlot(3).axisType = 'semilogy';
    postPlot(3).yname = 'error_p';

    postPlot(4) = postPlot(1);
    postPlot(5) = postPlot(2);
    postPlot(6) = postPlot(3);
    for j = 4:6
        postPlot(j).xname = 'f_ROM';
    end
    switch misc.BC
        case {'SHBC','NBC'}
            noDomains = 1;
        case 'SSBC'
            noDomains = 2;
        case 'NNBC'
            noDomains = 3;
    end
    varCol = setHetmaniukParameters(noDomains);
    msh.meshFile = 'createNURBSmesh_EL';
    msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
    msh.degree = 3;
    msh.M = 6; % 6
    if strcmp(misc.method{1},'PML')
        misc.formulation = {'GSB'};
        if msh.parm == 1
            varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-4)-1, max(2^(M-4)-1,3)];
%             varCol{1}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, 2^(M-1)/8-1, max(2^(M-4)-1,3)];
        else
            varCol{1}.refinement = @(M) [3*2^(M-3)-1, 3*2^(M-3)-1, 2^(M-1)/8-1, max(2^(M-4)-1,3)];
        end
    else
        misc.formulation = {'BGC'};
        varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-4)-1, max(2^(M-4)-1,3)];
    end
    misc.extraGP = [9-msh.degree,0,0];    % extra quadrature points
    if noDomains > 1
        if msh.parm == 1
            varCol{2}.refinement = @(M,t,t_fluid) [0, 2^(M-1)-1, max(round(t/t_fluid)*2^(M-1),0)];
        else
            varCol{2}.refinement = @(M,t,t_fluid) [3*2^(M-3)-1, 3*2^(M-3)-1, max(round(t/t_fluid)*2^(M-1),0)];
        end
    end
    if noDomains > 2
        varCol{3}.refinement = @(M) [0, 2^(M-1)-1, max(2^(M-1)-1,0)];
    end
    switch misc.BC
        case {'SHBC','NBC'}
            k = linspace(9, 36, 3);
%             k = linspace(9, 36, 3)/5;
            k_ROM = k(1):0.05:k(end);
%             k_ROM = k(1):0.2:k(end);
%             k_ROM = linspace(9, 36, 5);
            c_f = varCol{1}.c_f;
            rom.omega_ROM = k_ROM*c_f;
            f = k*c_f/(2*pi);
            omega = 2*pi*f;
            if hetmaniukCase
                misc.P_inc = -1;
                ffp.beta = pi/2;   
            end
            ffp.paramPts = {[0,1,0], []};
            postPlot(1).plotResults = true;
            postPlot(2).plotResults = true;
            postPlot(3).plotResults = true;
            postPlot(4).plotResults = false;
            postPlot(5).plotResults = false;
            postPlot(6).plotResults = false;
        case {'SSBC','NNBC'}
            f = linspace(1430, 4290, 5);
            omega = 2*pi*f;
            f_ROM = f(1):12:f(end);
%             f_ROM = f(1):120:f(end);
%             f_ROM = linspace(1430, 4290, 9);
            rom.omega_ROM = 2*pi*f_ROM;
            k = omega/varCol{1}.c_f;
            if hetmaniukCase
                misc.P_inc = 1;
                ffp.beta = -pi/2;   
            end
            ffp.paramPts = {[0,0,0], []};
            postPlot(1).plotResults = false;
            postPlot(2).plotResults = false;
            postPlot(3).plotResults = false;
            postPlot(4).plotResults = true;
            postPlot(5).plotResults = true;
            postPlot(6).plotResults = true;
    end
    misc.omega = omega;

    rom.basisROMcell = {'Pade','Taylor','DGP','Hermite','Bernstein'};  % do not put basisROMcell in loopParameters (this is done automatically)
    rom.basisROMcell = {'Pade','DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
%     rom.basisROMcell = {'DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
    iem.N = 16; % 9
    iem.p_ie = 5;
    iem.s_ie = 2;
    iem.IElocSup = 1;        % Toggle usage of radial shape functions in IE with local support
    

    pml.eps = 1e9*eps;      % choosing eps = eps yields machine precicion at Gamma_b, but requires more "radial" elements in the PML to resolve the rapid decay function
    pml.sigmaType = 1;  % sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n
    pml.dirichlet = false;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
    switch misc.BC
        case {'SHBC','NBC'}
            pml.t = 0.2;         % thickness of PML
            pml.gamma = 2.5;          % parameter for sigmaType = 1
            misc.r_a = 1.2;
            rom.noVecsArr = [8,16,24,32];
%             rom.noVecsArr = 32;
        case {'SSBC','NNBC'}
            rom.noVecsArr = [4,8,12,20];
%             rom.noVecsArr = 20;
            pml.t = 0.4;         % thickness of PML
            pml.gamma = 2.0;          % parameter for sigmaType = 1
            misc.r_a = 1.4;
    end
    pml.C = -log(pml.eps)*(pml.n+1)/(k(1)*pml.t);
    
    rom.useROM = true;

    misc.storeFullVarCol = false;
    misc.scatteringCase = 'Sweep';
    loopParameters = {'msh.M','misc.method','misc.coreMethod','misc.BC'};
    
    collectIntoTasks

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rom.useROM = false;
    for j = 1:6
        postPlot(j).noXLoopPrms   	= 1;
        postPlot(j).xname = postPlot(j).xname(1);
        postPlot(j).xLoopName     	= 'misc.omega';
    end
    misc.omega = rom.omega_ROM;
    
    misc.scatteringCase = 'BI';
    loopParameters = {'msh.M','misc.method','misc.coreMethod','misc.BC','misc.omega'};
    collectIntoTasks
    
    misc.omega = rom.omega_ROM(end);
    para.plotResultsInParaview = true;
    para.extraXiPts              = '1';  % Extra visualization points in the xi-direction per element
    para.extraEtaPts             = '1';  % Extra visualization points in the eta-direction per element
    para.extraZetaPts            = '1';   % Extra visualization points in the zeta-direction per element
%     collectIntoTasks

    misc.omega = rom.omega_ROM;
    para.plotResultsInParaview = 0;
    misc.method = {'BA'};
%     misc.formulation = {'SL2E'};
    misc.formulation = {'VL2E'};
    collectIntoTasks
end
