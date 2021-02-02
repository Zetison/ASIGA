function studies = getTask_Hetmaniuk()
% This study is based on Hetmaniuk2012raa (Figure 7)
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271

counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
hetmaniukCase = false; % evaluating solution at boundary not implemented

scatteringCase = 'Sweep';
% BC = 'NNBC';
model = 'S1';  % Spherical shell
coreMethod = {'IGA','hp_FEM'};
coreMethod = {'IGA'};
method = {'IE'};
% formulation = {'PGC'};
formulation = {'BGC'};
% BCs = {'SHBC'};
% BCs = {'SSBC'};
BCs = {'NNBC'};
% BCs = {'SHBC','SSBC'};
plotFarField = ~hetmaniukCase;

calculateFarFieldPattern    = true;     % Calculate far field pattern
alpha_s = 0;                            % Aspect angle of incident wave
beta_s  = -pi/2;                        % Elevation angle of incident wave
alpha   = 0;                            % Aspect angles of observation points
beta = -pi/2;   
r = 1;                            % radii for near-field evaluation.

r_a = 1.2;

parm = 1;
calculateSurfaceError = 0;
calculateVolumeError  = 1;
calculateFarFieldPattern = 1;
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;
% prePlot.colorFun = @(v) abs(norm2(v)-1);
prePlot.resolution = [20,20,0];
computeCondNumber = 0;

postPlot(1).yname        	= 'energyError';
% postPlot(1).yname        	= 'surfaceError';
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

for BC = BCs
    postPlot(1).xname = 'k_ROM';
    postPlot(2).xname = 'k_ROM';
    postPlot(3).xname = 'k_ROM';
    postPlot(4).xname = 'f_ROM';
    postPlot(5).xname = 'f_ROM';
    postPlot(6).xname = 'f_ROM';
    switch BC{1}
        case 'SHBC'
            noDomains = 1;
        case 'SSBC'
            noDomains = 2;
        case 'NNBC'
            noDomains = 3;
    end
    varCol = setHetmaniukParameters(noDomains);
    varCol{1}.meshFile = 'createNURBSmesh_EL';
    Xi = [0,0,0,1,1,2,2,3,3,3]/3;
    varCol{1}.refinement = @(M) [0, 2^(M-1)-1, max(2^(M-1)/8-1,0)];
    if noDomains > 1
        varCol{2}.refinement = @(M,t,t_fluid) [0, 2^(M-1)-1, max(round(t/t_fluid)*2^(M-1),0)];
    end
    if noDomains > 2
        varCol{3}.refinement = @(M) [0, 2^(M-1)-1, max(2^(M-1)-1,0)];
    end
    switch BC{1}
        case 'SHBC'
            k = linspace(9, 36, 3);
            k = linspace(9, 36, 3)/10;
            k_ROM = k(1):0.05:k(end);
%             k_ROM = k(1):0.5:k(end);
            c_f = varCol{1}.c_f;
            omega_ROM = k_ROM*c_f;
            f = k*c_f/(2*pi);
            if hetmaniukCase
                P_inc = -1;
                beta = pi/2;   
            end
            postPlot(1).plotResults = true;
            postPlot(2).plotResults = true;
            postPlot(3).plotResults = true;
            postPlot(4).plotResults = false;
            postPlot(5).plotResults = false;
            postPlot(6).plotResults = false;
        case {'SSBC','NNBC'}
            f = linspace(1430, 4290, 5);
            f = linspace(143, 429, 5);
            f_ROM = f(1):12:f(end);
%             f_ROM = f(1):120:f(end);
            omega_ROM = 2*pi*f_ROM;
            postPlot(1).xname = 'f_ROM';
            if hetmaniukCase
                P_inc = 1;
                beta = -pi/2;   
            end
            postPlot(1).plotResults = false;
            postPlot(2).plotResults = false;
            postPlot(3).plotResults = false;
            postPlot(4).plotResults = true;
            postPlot(5).plotResults = true;
            postPlot(6).plotResults = true;
    end
    omega = 2*pi*f;



    % noVecsArr = [4,8,16];        % do not put noVecsArr in loopParameters (this is done automatically)
    % basisROMcell = {'Bernstein','Hermite','Pade','Taylor','DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
    % basisROMcell = {'DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
    basisROMcell = {'Pade','Taylor','DGP','Hermite','Bernstein'};  % do not put basisROMcell in loopParameters (this is done automatically)
    % basisROMcell = {'Taylor'};  % do not put basisROMcell in loopParameters (this is done automatically)
    basisROMcell = {'DGP'};  % do not put basisROMcell in loopParameters (this is done automatically)
    % basisROMcell = {'Hermite'};  % do not put basisROMcell in loopParameters (this is done automatically)
    noVecsArr = 64;
    % noVecsArr = [2,4,8,16,32,64];
    % noVecsArr = 8;
    % k_start = 9/10;
    % k_end = 36/10;
    % P = 3;
    % n = P*noVecs;
    % temp = cos(pi/(2*n));
    % a = ((k_start+k_end)*temp-(k_end-k_start))/(2*temp);
    % b = ((k_start+k_end)*temp+(k_end-k_start))/(2*temp);
    % j = 1:n;
    % k_arr = 1/2*(a+b)+1/2*(b-a)*cos((2*n-2*j+1)*pi/(2*n));
    % k = k_arr(round(linspace(1,n,P)));
    degree = 4;
    % degree = 2:5;
    M = 3:5; % 5
    N = 7; % 9
    % M = 1; 
    % N = 2;
    useROM = true;
    % 
    % useROM = false;
    % postPlot(1).xname = 'k';
    % k = k_ROM;
    % c_f = 1500;
    % f = k*c_f/(2*pi);
    % refineThetaOnly = 0;
    p_ie = 5;
    s_ie = 2;
    IElocSup = 0;        % Toggle usage of radial shape functions in IE with local support

    storeFullVarCol = false;
    loopParameters = {'M','method','BC'};
    collectIntoTasks

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    useROM = false;
    postPlot(1).xname = postPlot(1).xname(1);
    postPlot(2).xname = postPlot(2).xname(1);
    postPlot(3).xname = postPlot(3).xname(1);
    postPlot(4).xname = postPlot(4).xname(1);
    postPlot(5).xname = postPlot(5).xname(1);
    postPlot(6).xname = postPlot(6).xname(1);
    omega = omega_ROM;
    f = omega/(2*pi);
    coreMethod = {'IGA'};
    method = {'BA'};
%     formulation = {'SL2E'};
    formulation = {'VL2E'};
    collectIntoTasks
    
    coreMethod = {'IGA'};
    method = {'IE'};
    formulation = {'BGU'};
    collectIntoTasks
end
