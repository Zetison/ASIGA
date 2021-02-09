function layer = extract_e3Dss_data(varCol)
M = numel(varCol);
layer = cell(1,M);
copyFields = {'P_inc','k','applyLoad','r_s','d_vec','isSphericalShell','N_max','omega','BC','a','noRHSs','model','splitExteriorFields'};

for i = 1:numel(copyFields)
    if isfield(varCol{1},copyFields{i})
        layer{1}.(copyFields{i}) = varCol{1}.(copyFields{i});
    end
end
for m = 1:M
    if isfield(varCol{m},'R_i')
        layer{m}.R_i = varCol{m}.R_i; % Inner radius of layer
    end
    layer{m}.rho              = varCol{m}.rho;    % Mass density
    layer{m}.media = varCol{m}.media;
    switch varCol{m}.media
        case 'fluid'
            layer{m}.c_f          	= varCol{m}.c_f;       % Speed of sound
            layer{m}.calc_p_0       = m == 1;      % Toggle calculation of the far field pattern
            layer{m}.calc_p       	= true;      % Toggle calculation of the scattered pressure
            layer{m}.calc_dp      	= true(1,3); % Toggle calculation of the three components of the gradient of the pressure
            layer{m}.calc_p_inc     = true;      % Toggle calculation of the incident pressure
            layer{m}.calc_dp_inc    = true(1,3); % Toggle calculation of the three components of the gradient of the incident pressure
        case {'solid','viscoelastic'}
            layer{m}.E            = varCol{m}.E;      % Youngs modulus for solid layers
            layer{m}.nu           = varCol{m}.nu;        % Poisson ratio for solid layers
            layer{m}.calc_u       = true(1,3); % Toggle calculation of the three components of the displacement
            layer{m}.calc_du      = true(3,3); % Toggle calculation of the three cartesian derivatives of the three components of the displacement [du_xdx du_xdy du_xdz; 
                                                %                                                                                                    du_ydx du_ydy du_ydz; 
                                                %                                                                                                    du_zdx du_zdy du_zdz]
            layer{m}.calc_sigma   = true(1,6); % Toggle calculation of the six components of the stress field (cartesian coordinates) [sigma_xx sigma_yy sigma_zz sigma_yz sigma_xz sigma_xy]
    end
end