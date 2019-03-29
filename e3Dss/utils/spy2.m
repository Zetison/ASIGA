function spy2(A)

B = log10(abs(full(A)));
B(abs(full(A)) == 0) = Inf;
if any(isinf(B(:)))
    B(B == Inf) = min(min(B))+min(min(B))/5;
    issparse = true;
else
    issparse = false;
end
% infRGB = reshape([1 1 1],1,1,3);   
% infimgdata = repmat(infRGB, size(B,1), size(B,2)); 
% image(infimgdata, 'alphadata', ~isnan(B));  
% hold on
% imagesc(B,'alphadata', ~(isnan(B)|isinf(B)));
imagesc(B);
% colormap(jet());
% map = [linspace(0,1,1000).', zeros(1000,1), linspace(1,0,1000).'];
% colormap(map)
if 1
    map = colormap;
    if issparse
        map = [1,1,1; map];
    end
    colormap(map)
else
    % colormap gray
    map = colormap;
    map = flipud(map);
    colormap(map)
end
% colormap default
% colormap autumn
colorbar
hold off