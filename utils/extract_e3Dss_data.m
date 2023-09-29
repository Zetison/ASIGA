function layer = extract_e3Dss_data(task)
M = numel(task.varCol);
layer = cell(1,M);
copyFields = {'P_inc','applyLoad','r_s','d_vec','isSphericalShell','N_max','omega','BC','a','noRHSs','model','splitExteriorFields'};

for i = 1:numel(copyFields)
    if isfield(task.misc,copyFields{i})
        layer{1}.(copyFields{i}) = task.misc.(copyFields{i});
    elseif isfield(task,copyFields{i})
        layer{1}.(copyFields{i}) = task.(copyFields{i});
    elseif isfield(task.varCol{1},copyFields{i})
        layer{1}.(copyFields{i}) = task.varCol{1}.(copyFields{i});
    end
end
for m = 1:M
    if isfield(task.varCol{m},'R')
        layer{m}.R = task.varCol{m}.R; % Inner radius of layer
    end
    layer{m}.rho   = task.varCol{m}.rho;    % Mass density
    layer{m}.media = task.varCol{m}.media;
    switch task.varCol{m}.media
        case 'fluid'
            layer{m}.c_f          	= task.varCol{m}.c_f;       % Speed of sound
            layer{m}.calc_p_0       = m == 1;      % Toggle calculation of the far field pattern
            layer{m}.calc_p       	= true;      % Toggle calculation of the scattered pressure
            layer{m}.calc_dp      	= true(1,3); % Toggle calculation of the three components of the gradient of the pressure
            layer{m}.calc_p_inc     = true;      % Toggle calculation of the incident pressure
            layer{m}.calc_dp_inc    = true(1,3); % Toggle calculation of the three components of the gradient of the incident pressure
        case {'solid','viscoelastic'}
            layer{m}.E            = task.varCol{m}.E;      % Youngs modulus for solid layers
            layer{m}.nu           = task.varCol{m}.nu;        % Poisson ratio for solid layers
            layer{m}.calc_u       = true(1,3); % Toggle calculation of the three components of the displacement
            layer{m}.calc_du      = true(3,3); % Toggle calculation of the three cartesian derivatives of the three components of the displacement [du_xdx du_xdy du_xdz; 
                                                %                                                                                                    du_ydx du_ydy du_ydz; 
                                                %                                                                                                    du_zdx du_zdy du_zdz]
            layer{m}.calc_sigma   = true(1,6); % Toggle calculation of the six components of the stress field (cartesian coordinates) [sigma_xx sigma_yy sigma_zz sigma_yz sigma_xz sigma_xy]
    end
end