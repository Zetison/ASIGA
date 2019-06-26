scatteringCase = 'Sweep'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'IL';
formulation = 'BGU';
method = 'IE';

coreMethods = {'IGA','IGA','hp_FEM','linear_FEM'}; % [5, 4, 2, 1, 3]
coreMethods = {'IGA'}; % [5, 4, 2, 1, 3]
for i_coreM = 1:length(coreMethods) %{'IGA'}
    coreMethod = {coreMethods{i_coreM}};
    for BC = {'SSBC'} %{'SHBC','SSBC','NNBC'} %
        k = 1;
        c_f = 1524;
        omega = k*c_f;
        f = omega/(2*pi); 

        switch coreMethod{1}
            case 'IGA'
                if i_coreM == 1
                    M = 4; % 4
                    degreeElev = 1; %1
                else
                    M = 5;
                    degreeElev = 0:1;
                end
            case 'hp_FEM'
                degreeElev = 0;
                M = 5;
            case 'linear_FEM'
                degreeElev = 0;
                M = 6;
        end

        alpha_s = pi;
        beta_s = 0;
        alpha = alpha_s;
        beta = beta_s;

        plotResultsInParaview = 0;
        calculateFarFieldPattern = 1;
        calculateVolumeError = 1;
        calculateSurfaceError = 0;
        loopParameters = {'M','degree','coreMethod','BC'};

        collectIntoTasks
    end
end