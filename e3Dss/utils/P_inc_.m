function p = P_inc_(omega,omega_c,type)

switch type
    case 1
        p = 4/sqrt(3)*omega_c^3./((omega.^2-omega_c^2).*(omega.^2-4*omega_c^2)).*(1-exp(-2*pi*1i*omega/omega_c));
        indices = logical( (abs(abs(omega)-abs(omega_c))/abs(omega_c) < 10*eps) ...
                          +(abs(abs(omega)-abs(2*omega_c))/abs(2*omega_c) < 10*eps));
        p(indices) = 4/(3*sqrt(3))*1i*pi./omega(indices).*exp(1i*pi*omega(indices)/omega_c);
    case 2
        p = -3/2*1i*omega*omega_c^2./((omega.^2-omega_c^2).*(omega.^2-4*omega_c^2)).*(1-exp(-2*pi*1i*omega/omega_c));
        indices = logical( (abs(abs(omega)-abs(omega_c))/abs(omega_c) < 10*eps) ...
                          +(abs(abs(omega)-abs(2*omega_c))/abs(2*omega_c) < 10*eps));
        p(indices) = pi/(2*omega_c)*exp(1i*pi*omega(indices)/omega_c);
end