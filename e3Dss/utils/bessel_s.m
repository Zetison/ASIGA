function Z = bessel_s(n,z,type)
%Returns the n'th spherical bessel function of kind "type" evaluated at
%every element in z

if isa(z,'sym')
    PI = vpa('pi');
    tiny = realmin('double');
elseif isa(z,'mp')
    PI = mp('pi');
    tiny = realmin('double');
else
    tiny = eps;
    PI = pi;
end
if type == 1
    indices = logical(abs(z) < tiny);
    z(indices) = 1; % Avoid symbolic division by zero. These values will be replaced anyways
end

Z = sqrt(PI/2)./sqrt(z).*bessel_c(n+0.5,z,type);

if type == 1
    if n == 0
        Z(indices) = 1;
    else
        Z(indices) = 0;
    end
end
if ~(isa(z,'sym') || isa(z,'mp'))
    if any(any(abs(Z) > 1e290))
        error('A Bessel function evaluation was too large')
    end
end