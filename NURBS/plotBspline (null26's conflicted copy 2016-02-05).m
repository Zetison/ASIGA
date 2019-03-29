function h = plotBspline(i,p,n,Xi,noPts,plotDers)

if nargin < 6
    % As findKnotSpan will give a new set of functions at xi = Xi(i+p+1) we
    % skip this last point, and add it later
    xi_array = linspace2(Xi(i), Xi(i+p+1), noPts);

    B = zeros(1,length(xi_array));

    for c = 1:noPts
        xi = xi_array(c);
        i1 = findKnotSpan(n, p, xi, Xi);
        N = Bspline_basis(i1, xi, p, Xi);
        B(c) = N(p+1+i-i1);
    end
    if Xi(1) == Xi(i+p)
        xi_array = [Xi(1) xi_array Xi(i+p+1) Xi(end)];
        B = [1 B 0 0];
    elseif Xi(i+1) == Xi(end)
        xi_array = [Xi(1) Xi(i) xi_array Xi(end)];
        B = [0 0 B 1];
    else
        xi_array = [Xi(1) Xi(i) xi_array Xi(i+p+1) Xi(end)];
        B = [0 0 B 0 0];
    end
    h = plot(xi_array, B);
else
    % As findKnotSpan will give a new set of functions at xi = Xi(i+p+1) we
    % skip this last point, and add it later
    xi_array = linspace2(Xi(i), Xi(i+p+1), noPts);

    B = zeros(1,length(xi_array));
    Bders = zeros(1,length(xi_array));

    for c = 1:noPts
        xi = xi_array(c);
        i1 = findKnotSpan(n, p, xi, Xi);
        [N, dNdxi] = Bspline_basisDers(i1, xi, p, Xi);
        B(c) = N(p+1+i-i1);
        Bders(c) = dNdxi(p+1+i-i1);
    end
    if Xi(1) == Xi(i+p)
        xi_array = [Xi(1) xi_array Xi(i+p+1) Xi(end)];
        B = [1 B 0 0];
        Bders = [0 Bders 0 0];
    elseif Xi(i+1) == Xi(end)
        xi_array = [Xi(1) Xi(i) xi_array Xi(end)];
        B = [0 0 B 1];
        Bders = [0 0 Bders 0];
    else
        xi_array = [Xi(1) Xi(i) xi_array Xi(i+p+1) Xi(end)];
        B = [0 0 B 0 0];
        Bders = [0 0 Bders 0 0];
    end
    h = plot(xi_array, B, xi_array, Bders);
end

