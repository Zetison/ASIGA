M = length(R_o);
if SHBC
    E = E(1:end-1);
    rho_s = rho_s(1:end-1);
    rho_f = rho_f(1:end-1);
    c_f = c_f(1:end-1);
    nu = nu(1:end-1);
    R_i = R_i(1:end-1);
    BC = 'SHBC';
elseif ESBC
    rho_f = rho_f(1:end-1);
    c_f = c_f(1:end-1);
    R_i = R_i(1:end-1);
    BC = 'ESBC';
elseif SSBC
    rho_f = rho_f(1:end-1);
    c_f = c_f(1:end-1);
    BC = 'SSBC';
else
    BC = 'NNBC';
end