function studies = getTask_articleIGA_IL_sweep()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 12-17 in Venas2018iao
% Venas2018iao is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
% misc.scatteringCase = 'Sweep'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'IL';
formulation = 'BGU';
misc.method = 'IE';

axisymmetric = 1;
if axisymmetric
    extraGP = [7,0,0];    % extra quadrature points
    alpha_s = 0;
    beta_s = -pi/2;
    Xi = [0,0,0,1,1,2,2,3,3,3]/3;
    refineThetaOnly = true;
else
    alpha_s = pi;
    beta_s = 0;
end
alpha = alpha_s;
beta = beta_s;
parm = 1;
N = 4; % not sure if this was the original choice in the article

prePlot.plot2Dgeometry = 0;  % Plot cross section of mesh and geometry
prePlot.plot3Dgeometry = 0;  % Plot visualization of mesh and geometry in 3D
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting

runTasksInParallel = 0;

postPlot(1).xname        	= 'k';
postPlot(1).yname        	= 'TS';
postPlot(1).lineStyle    	= '-';
if runTasksInParallel || strcmp(misc.scatteringCase,'BI')
    postPlot(1).xLoopName     	= 'f';
    postPlot(1).noXLoopPrms   	= 1;
end

postPlot(2) = postPlot(1);
postPlot(2).axisType      	= 'semilogy';
postPlot(2).yname        	= 'energyError';

plotResultsInParaview = 0;
calculateFarFieldPattern = 1;
calculateVolumeError = 1;
err.calculateSurfaceError = 0;
% extraGP = -2;
if runTasksInParallel
    loopParameters = {'M','degree','misc.coreMethod','BC','f'};
else
%     loopParameters = {'M','degree','misc.coreMethod','BC'};
    loopParameters = {'M','degree','misc.coreMethod','BC','f'};
end
BCs = {'SHBC','SSBC','NNBC'};
% BCs = {'NNBC'};
% BCs = {'SHBC'};
misc.coreMethods = {'IGA','IGA','linear_FEM','hp_FEM'}; % [5, 4, 2, 1, 3]
% misc.coreMethods = {'linear_FEM'}; % [5, 4, 2, 1, 3]
M_0 = 4;
% misc.coreMethods = {'IGA'}; % [5, 4, 2, 1, 3]
for i_coreM = 1:length(misc.coreMethods) %{'IGA'}
    misc.coreMethod = {misc.coreMethods{i_coreM}};
    for BC = BCs %
        npts = 1500;
%         npts = 15;
%         npts = 2;
        if strcmp(BC{1},'SSBC')
            npts = npts*2;
            specialValues = [0.317054564603519 %
                               0.392367003924396 %
                               0.450625070843336 %
                               0.495802583768421 %
                               0.538194539627076 %
                               0.584225864137528 %
                               0.638340507928227 %
                               0.703448659735311 %
                               0.781278648263076 %
                               0.872642563605135 %
                               0.977700430550731 %
                               1.096200471571476 %
                               1.227664074628201 %
                               1.371508743070356 %
                               1.527120122444878 %
                               1.693889192704899 %
                               1.871228417747961 %
                               ];
            noDomains = 2;
        elseif strcmp(BC{1},'NNBC')
            npts = npts*2;
            specialValues = [0.250621182794531 %
                           0.320300579445871 %
                           0.370671479527136 %
                           0.412992731010227 %
                           0.454191270410376 %
                           0.499088778976889 %
                           0.551286412239756 %
                           0.613370456080303 %
                           0.687008309336546 %
                           0.773084257718564 %
                           0.871890313908958 %
                           0.983323027396819 %
                           1.107045032953710 %
                           1.242597693362253 %
                           1.389470517759271 %
                           1.547139666101034 %
                           1.715087015757087 %
                           1.892808062465205 %
                           ];
            noDomains = 3;
        else
            specialValues = [];
            noDomains = 1;
        end
        varCol = setIhlenburgParameters(noDomains);
        varCol{1}.meshFile = 'createNURBSmesh_EL';
        if axisymmetric
            Xi = [0,0,0,1,1,2,2,3,3,3]/3;
            varCol{1}.refinement = @(M) [0, 2^(M-1)-1, max(2^(M-1)/8-1,0)];
            if noDomains > 1
                varCol{2}.refinement = @(M,t,t_fluid) [0, 2^(M-1)-1, max(round(t/t_fluid)*2^(M-1),0)];
            end
            if noDomains > 2
                varCol{3}.refinement = @(M) [0, 2^(M-1)-1, max(2^(M-2)-1,0)];
            end
        else
            varCol{1}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, max(2^(M-1)/8-1,0)];
            if numel(varCol) > 1
                varCol{2}.refinement = @(M,t,t_fluid) [2^(M-1)-1, 2^(M-1)-1, max(round(t/t_fluid)*2^(M-1),0)];
            end
            if numel(varCol) > 2
                varCol{3}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, max(2^(M-2)-1,0)];
            end
        end
        postPlot(1).xScale = varCol{1}.R;
        postPlot(2).xScale = varCol{1}.R;
%         specialValues = [];
        k = linspace(0.001, 2, npts);
        k = unique(sort([k'; specialValues; 1]'));
%         k = [0.99, 1];
        omega = k*varCol{1}.c_f;
        f = omega/(2*pi); 

        switch misc.coreMethod{1}
            case 'IGA'
                if i_coreM == 1
                    M = (M_0-1):M_0; % 4
                    degree = 3; %3
                else
                    M = M_0+1;
                    degree = 2:3;
                end
            case 'hp_FEM'
                degree = 2;
                M = M_0+1;
            case 'linear_FEM'
                degree = 2;
                M = M_0+2;
        end
        collectIntoTasks
    end
end