%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 12-17 in Venas2018iao
% Venas2018iao is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'IL';
formulation = 'BGU';
method = 'IE';

alpha_s = pi;
beta_s = 0;
alpha = alpha_s;
beta = beta_s;
parm = 1;
N = 4; % not sure if this was the original choice in the article

prePlot.plot2Dgeometry = 0;  % Plot cross section of mesh and geometry
prePlot.plot3Dgeometry = 0;  % Plot visualization of mesh and geometry in 3D
plotResultsInParaview = 0;
calculateFarFieldPattern = 1;
calculateVolumeError = 1;
calculateSurfaceError = 0;
runTasksInParallel = 0;
extraGP = -2;
loopParameters = {'M','degree','coreMethod','BC','f'};

% for coreMethod = {'IGA','hp_FEM','linear_FEM'}
%     for BC = {'SHBC','SSBC','NNBC'}
coreMethods = {'linear_FEM','hp_FEM','IGA','IGA'}; % [5, 4, 2, 1, 3]
% coreMethods = {'linear_FEM'}; % [5, 4, 2, 1, 3]
coreMethods = {'IGA'}; % [5, 4, 2, 1, 3]
for i_coreM = 1:length(coreMethods) %{'IGA'}
    coreMethod = {coreMethods{i_coreM}};
    for BC = {'NNBC','SSBC','SHBC'} %
%         npts = 1500;
        npts = 2;
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
%         specialValues = [];
        k = linspace(0.001, 2, npts);
        c_f = 1524;
        k = unique(sort([k'; specialValues; 1]'));
%         k = [0.99, 1];
        omega = k*c_f;
        f = omega/(2*pi); 

        switch coreMethod{1}
            case 'IGA'
                if i_coreM == 3
                    M = 4; % 4
                    degree = 3; %1
                else
                    M = 5; % 5
                    degree = 2;
                end
            case 'hp_FEM'
                degree = 2;
                M = 5; % 5
            case 'linear_FEM'
                degree = 1;
                M = 6; % 6
        end
        collectIntoTasks
    end
end