function d2Z = d2bessel_s(n,z,i,Z)
% Returns the second derivative of the n'th spherical bessel function of kind i
% evaluated at every element in z

if isa(z,'sym') || isa(z,'mp')
    tiny = realmin('double');
else
    tiny = eps;
end
if i == 1
    indices = logical(abs(z) < tiny);
    z(indices) = 1; % Avoid symbolic division by zero. These values will be replaced anyways
end
d2Z = (n*(n-1)./z.^2-1).*Z{i,1} + 2./z.*Z{i,2};
if i == 1
    if n == 0
        d2Z(indices) = -ones(1,class(z))/3;
    elseif n == 2
        d2Z(indices) = 2*ones(1,class(z))/15;
    else
        d2Z(indices) = 0;
    end
end
if ~(isa(z,'sym') || isa(z,'mp'))
    if any(any(abs(d2Z) > 1e290))
        error('A Bessel function evaluation was too large')
    end
end