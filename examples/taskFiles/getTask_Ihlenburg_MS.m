function studies = getTask_Ihlenburg_MS()
% This study is based on Ihlenburg1998fea (Figure 7)
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271
% Experimental data at scale 1:50
counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
scatteringCase = 'Sweep';
BC = 'SSBC';
model = 'IMS';  % Spherical shell

coreMethod = {'IGA','hp_FEM'};
coreMethod = {'IGA'};
method = {'IE'};
applyLoad = 'pointSource';
applyLoad = 'planeWave';
formulation = {'BGC'};
BCs = {'SSBC'};
% BCs = {'SHBC','SSBC'};
% BCs = {'SHBC'};

warning('off','NURBS:weights')

alpha_s = 0;
beta_s = 0;   
alpha   = alpha_s;                            % Aspect angles of observation points
beta = beta_s;   

prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;
% prePlot.colorFun = @(v) abs(norm2(v)-1);
prePlot.resolution = [20,20,0];

postPlot(1).yname        	= 'TS';
% postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'plot';
postPlot(1).lineStyle    	= '-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).legendEntries 	= {};
postPlot(1).subFolderName 	= '';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;

for BC = BCs
    postPlot(1).xname = 'k_ROM';
    switch BC{1}
        case 'SHBC'
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

    varCol{1}.meshFile = 'createNURBSmesh_M3';
    k = linspace(2.5, 20, 5)/varCol{1}.R1;

    k_ROM = k(1):0.05:k(end);
    c_f = varCol{1}.c_f;
    omega_ROM = k_ROM*c_f;
    f = k*c_f/(2*pi);
    omega = 2*pi*f;
    
    Xi = [0,0,0,1,1,2,2,3,3,3]/3;
    varCol{1}.refinement = @(M) [0, 2^(M-1)-1, max(2^(M-1)/8-1,0)];
    if noDomains > 1
        varCol{2}.refinement = @(M,t,t_fluid) [0, 2^(M-1)-1, max(round(t/t_fluid)*2^(M-1),0)];
    end
    
    sdfmsc = 3*(varCol{1}.R1+varCol{1}.L+varCol{1}.R2)/1.36;
    r_s = sdfmsc - varCol{1}.L/2;
    c_z = 45;
    c_xy = 15;

    basisROMcell = {'Pade','Taylor','DGP','Hermite','Bernstein'};  % do not put basisROMcell in loopParameters (this is done automatically)
    basisROMcell = {'DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
    noVecsArr = 64;
    degree = 2;
    M = 1; % 5
    N = 7; % 9
    useROM = true;
    p_ie = 5;
    s_ie = 2;
    IElocSup = 0;        % Toggle usage of radial shape functions in IE with local support

    storeFullVarCol = false;
    loopParameters = {'M','method','BC'};
    collectIntoTasks

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    useROM = false;
    postPlot(1).xname = postPlot(1).xname(1);
    omega = omega_ROM;
    f = omega/(2*pi);
    coreMethod = {'IGA'};
    method = {'BA'};
%     formulation = {'SL2E'};
    formulation = {'VL2E'};
%     collectIntoTasks
    
    coreMethod = {'IGA'};
    method = {'IE'};
    formulation = {'BGU'};
%     collectIntoTasks
end
