function dZ = dbessel_s(n,z,i,Z,scaled)
% Returns the derivative of the n'th spherical bessel function of kind i
% evaluated at every element in z
if nargin < 5
    scaled = false;
end
if scaled
    if iscell(Z)
        dZ = n*Z{i,1} - z.*Z{i,2};
    else
        dZ = n*bessel_s(n,z,i) - z.*bessel_s(n+1,z,i);
    end
else
    if iscell(Z)
        if isa(z,'sym') || isa(z,'mp')
            tiny = realmin('double');
        else
            tiny = eps;
        end
        if i == 1
            indices = logical(abs(z) < tiny);
            z(indices) = 1; % Avoid symbolic division by zero. These values will be replaced anyways
        end
        dZ = n./z.*Z{i,1} - Z{i,2};
        if i == 1
            if n == 1
                dZ(indices) = ones(1,class(z))/3;
            else
                dZ(indices) = 0;
            end
        end
    else
        if isa(z,'sym') || isa(z,'mp')
            tiny = realmin('double');
        else
            tiny = eps;
        end
        if i == 1
            indices = logical(abs(z) < tiny);
            z(indices) = 1; % Avoid symbolic division by zero. These values will be replaced anyways
        end
        dZ = n./z.*bessel_s(n,z,i) - bessel_s(n+1,z,i);
        if i == 1
            if n == 1
                dZ(indices) = ones(1,class(z))/3;
            else
                dZ(indices) = 0;
            end
        end
    end
end
if ~isa(z,'sym')
    if any(any(abs(dZ) > 1e290))
        error('A Bessel function evaluation was too large')
    end
end