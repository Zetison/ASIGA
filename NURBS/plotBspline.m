function h = plotBspline(i,p,n,Xi,noPts,plotDers)
if nargin < 6
    plotDers = false;
end
if ~plotDers
    % As findKnotSpan will give a new set of functions at xi = Xi(i+p+1) we
    % skip this last point, and add it later
    
    xi_array = sort(unique([Xi(and(Xi(i) < Xi,Xi < Xi(i+p+1))), linspace2(Xi(i), Xi(i+p+1), noPts)]));
%     if Xi(i+1) == Xi(i+p) && Xi(i+1) ~= Xi(end) && Xi(i+1) ~= Xi(1)
%         xi_array = sort([Xi(i+1) xi_array]);
%     end
    B = zeros(1,length(xi_array));

    for c = 1:length(xi_array)
        xi = xi_array(c);
        i1 = findKnotSpan(n, p, xi, Xi);
        N = Bspline_basis(i1, xi, p, Xi, 0);
        B(c) = N(p+1+i-i1);
    end
    if Xi(1) == Xi(i+p)
        xi_array = [Xi(1) xi_array Xi(i+p+1) Xi(end)];
        B = [1 B 0 0];
    elseif Xi(i+1) == Xi(end)
        xi_array = [Xi(1) Xi(i) xi_array Xi(end)];
        B = [0 0 B 1];
    elseif Xi(i) == Xi(i+p)        
        xi_array = [Xi(1) Xi(i) Xi(i) xi_array Xi(i+p+1) Xi(end)];
        B = [0 0 1 B 0 0];
    elseif Xi(i+1) == Xi(i+p+1)        
        xi_array = [Xi(1) Xi(i) xi_array Xi(i+p+1) Xi(i+p+1) Xi(end)];
        B = [0 0 B 1 0 0];
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
    Bders2 = zeros(1,length(xi_array));
    Bders3 = zeros(1,length(xi_array));

    for c = 1:noPts
        xi = xi_array(c);
        i1 = findKnotSpan(n, p, xi, Xi);
        [N, dNdxi, d2Ndxi2, d3Ndxi3] = Bspline_basisDers2(i1, xi, p, Xi);
        B(c) = N(p+1+i-i1);
        Bders(c) = dNdxi(p+1+i-i1);
        Bders2(c) = d2Ndxi2(p+1+i-i1);
        Bders3(c) = d3Ndxi3(p+1+i-i1);
    end
    if Xi(1) == Xi(i+p)
        xi_array = [Xi(1) xi_array Xi(i+p+1) Xi(end)];
        B = [1 B 0 0];
        Bders = [0 Bders 0 0];
        Bders2 = [0 Bders2 0 0];
        Bders3 = [0 Bders3 0 0];
    elseif Xi(i+1) == Xi(end)
        xi_array = [Xi(1) Xi(i) xi_array Xi(end)];
        B = [0 0 B 1];
        Bders = [0 0 Bders 0];
        Bders2 = [0 0 Bders2 0];
        Bders3 = [0 0 Bders3 0];
    else
        xi_array = [Xi(1) Xi(i) xi_array Xi(i+p+1) Xi(end)];
        B = [0 0 B 0 0];
        Bders = [0 0 Bders 0 0];
        Bders2 = [0 0 Bders2 0 0];
        Bders3 = [0 0 Bders3 0 0];
    end
%     h = plot(xi_array, B, xi_array, Bders/30, xi_array, Bders2/400, xi_array, Bders3/4000);
%     legend('N','dNdxi','d2Ndxi2','d3Ndxi3')
    h = plot(xi_array, B, xi_array, Bders, xi_array, Bders2, xi_array, Bders3);
    legend('N','dNdxi','d2Ndxi2','d3Ndxi3')
%     h = plot(xi_array, Bders2, xi_array, Bders3);
%     legend('d2Ndxi2','d3Ndxi3')
end

