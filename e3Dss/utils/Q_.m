function Q = Q_(j,theta,P,dP,d2P,scaleByCsc)
if nargin < 6
    scaleByCsc = false;
end
if ~mod(j,2) && scaleByCsc && (any(theta == 0) || any(theta == pi))
    error('Q is not defined in this case (when theta = 0 or theta = pi)')
end
switch j
    case 0
        Q = P;
    case 1
        if scaleByCsc
            Q = -dP;
        else
            Q = -sin(theta).*dP;
        end
    case 2
        Q = -cos(theta).*dP + sin(theta).^2.*d2P;
    otherwise
        error('Not implemented')
end



