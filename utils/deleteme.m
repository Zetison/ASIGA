function r = deleteme(zeta,r_a,r_b,type)

switch type
    case 1
        r = zeros(size(zeta));
        indices = zeta < 1;
        r(indices) = 1./((1/r_b-1/r_a)*zeta(indices)+1/r_a);
        indices = zeta > 1;
        r(indices) = r_b./zeta(indices);
    case 2
        r = r_a./(zeta+1);
    case 3
        r = zeros(size(zeta));
        indices = zeta < 1;
        r(indices) = 1./(0.5*(1/r_b-1/r_a)*zeta(indices)+0.5*(1/r_b+1/r_a));
        indices = zeta >= 1;
        r(indices) = 2*r_b./(1-zeta(indices));
end

