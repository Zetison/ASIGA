function studies = getTask_PML_convergenceAnalysis()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 9 and 10 in Venas2018iao
% Venas2018iao is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
misc.model = 'IL';
noCoresToUse = 12;

msh.meshFile = 'createNURBSmesh_EL';
msh.parm = 1;
msh.explodeNURBS = false;   % Create patches from all C^0 interfaces
misc.checkNURBSweightsCompatibility = false;
misc.computeCondNumber = 0;       % Compute the condition number of the global matrix
prePlot.plotGeometryInfo    = 1;       % Plot domain boundaries (i.e. Gamma, Gamma_a, Neumann, Dirichlet, ...)
% misc.method = 'BA';
% misc.method = {'IENSG'};
% misc.method = {'BEM'};
% BC = {'SHBC', 'SSBC','NNBC'};
% for BC = {'SHBC', 'SSBC','NNBC'}
axiSymmetricCase = 1;
if axiSymmetricCase
    ffp.alpha_s = 0;
    ffp.beta_s = -pi/2;
else
    ffp.alpha_s = pi;
    ffp.beta_s = 0;
end

prePlot.plot3Dgeometry   = 1;
prePlot.plotFullDomain   = 0;        % Plot volumetric domains
prePlot.view             = [0,0];
prePlot.plotSubsets      = {'xz'};
prePlot.plotControlPolygon  = 0;       % Plot the control polygon for the NURBS mesh
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
% prePlot.colorFun = @(v) abs(norm2(v)-(r_a+t_PML));
% prePlot.resolution = [100,100,0];
warning('off','NURBS:weights')

% postPlot(1).xname       	= 'nepw';
postPlot(1).xname       	= 'dofs';
postPlot(1).yname        	= 'energyError';
postPlot(1).plotResults  	= 1;
postPlot(1).printResults 	= 1;
postPlot(1).axisType        = 'loglog';
postPlot(1).lineStyle   	= '*-';
postPlot(1).xLoopName     	= 'msh.M';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).legendEntries 	= {'msh.degree','pml.sigmaType','pml.n','misc.method','misc.coreMethod','misc.formulation','pml.t'};
% postPlot(2) = postPlot(1);
% postPlot(2).yname       	= 'cond_number';
% msh.explodeNURBS = prePlot.plot3Dgeometry;

msh.meshFile = 'createNURBSmesh_EL';
msh.refineThetaOnly = ffp.beta_s == -pi/2;
    
if msh.refineThetaOnly 
    msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
else
    msh.Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
end
connectedParameters = {{'msh.M','iem.N'},{'pml.sigmaType','pml.n'}};
% misc.extraGP = [3,3,3];

for method = {'PML'} %,'IE','BA'}
    misc.method = method{1};
    switch method{1}
        case 'IE'
            misc.formulation = {'BGU'};
        case 'PML'
            misc.formulation = {'GSB'};
        case 'BA'
            misc.formulation = {'VL2E'};
    end

    M_max = 6; % 7
    for BC = {'SHBC'}
        misc.BC = BC{1};
        c_f = 1524;
        k = 1;
        misc.omega = k*c_f;
        msh.degree = 2;
        if strcmp(BC{1}, 'SHBC')
            varCol = setIhlenburgParameters(1);
            msh.M = 1:M_max; %1:7
        elseif strcmp(BC{1}, 'SSBC')
            varCol = setIhlenburgParameters(2);
            msh.M = 1:M_max-1; %1:6
        elseif strcmp(BC{1}, 'NNBC')
            varCol = setIhlenburgParameters(3);
            msh.M = 1:M_max-2; %1:5
        end
%         msh.M = (M_max-1):M_max; %1:5
        msh.M = 6;
        pml.dirichlet = true;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
        misc.r_a = 1.25*varCol{1}.R;
%         misc.r_a = 1.5*varCol{1}.R;
%         varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-4)-1, max(iem.N - msh.degree,2^(M-3)-1)];
%         varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-4)-1, 2^(M-4)-1];
        pml.refinement = @(M,t) round(t/(0.25*varCol{1}.R)*2^(M-4)-1);
        if msh.refineThetaOnly
            varCol{1}.refinement = @(M,t) [0, 2^(M-1)-1, 2^(M-4)-1];
        else
            varCol{1}.refinement = @(M,t) [2^(M-1)-1, 2^(M-1)-1, 2^(M-4)-1];
        end

        para.plotResultsInParaview = 0;
        ffp.calculateFarFieldPattern = 0;
        err.calculateVolumeError = 1;
        err.calculateSurfaceError = 1;
        loopParameters = {'msh.M','msh.degree','pml.sigmaType','misc.method','misc.coreMethod','misc.formulation','misc.BC','pml.t','misc.r_a'};


%         for coreMethod = {'IGA'}
        for coreMethod = {'IGA','hp_FEM','h_FEM','C0_IGA'}
            misc.coreMethod = coreMethod{1};
            if strcmp(coreMethod{1},'IGA')
                if strcmp(method{1}, 'PML')
                    pml.t = [0.5*varCol{1}.R, 0.25*varCol{1}.R];         % thickness of PML
                else
                    pml.t = 0.25*varCol{1}.R;         % thickness of PML
                end
                if strcmp(method{1},'PML')
                    pml.sigmaType = [3,5];   	% sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n, sigmaType = 3: sigma(xi) = C/(1-xi)^n
                    pml.n = [1,2];
                else
                    pml.sigmaType = 3;
                    pml.n = 1;
                end
            	iem.N = floor(abs(2.^(msh.M-4)-1)) + msh.degree;
            else
                pml.t = 0.25*varCol{1}.R;         % thickness of PML
                pml.sigmaType = 3;
                pml.n = 1;
            	iem.N = (floor(abs(2.^(msh.M-4)-1)) + 1)*msh.degree;
            end
%             collectIntoTasks
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for degree = 4 %3:4
            iem.N = floor(abs(2.^(msh.M-4)-1)) + degree;
            msh.degree = degree;
            misc.coreMethod = 'IGA';
            pml.t = 0.25*varCol{1}.R;         % thickness of PML
%             pml.t = [0.5*varCol{1}.R, 0.25*varCol{1}.R];         % thickness of PML
%             if strcmp(method{1}, 'PML')
%                 pml.t = [0.5*varCol{1}.R, 0.25*varCol{1}.R];         % thickness of PML
%             else
%                 pml.t = 0.25*varCol{1}.R;         % thickness of PML
%             end
            collectIntoTasks
        end 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        misc.coreMethod = 'linear_FEM';
        msh.degree = 1;
        iem.N = floor(abs(2.^(msh.M-4)-1)) + 1;

%         collectIntoTasks 
    end
end