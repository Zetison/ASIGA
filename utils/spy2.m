function spy2(A,blackAndWhite)

if nargin < 2
    blackAndWhite = false;
end
B = log10(abs(full(A)));
B(abs(full(A)) == 0) = Inf;
B(B == Inf) = min(min(B))+min(min(B))/5;
% infRGB = reshape([1 1 1],1,1,3);   
% infimgdata = repmat(infRGB, size(B,1), size(B,2)); 
% image(infimgdata, 'alphadata', ~isnan(B));  
% hold on
% imagesc(B,'alphadata', ~(isnan(B)|isinf(B)));
imagesc(B);
axis image
% colormap(jet());
if blackAndWhite
    colormap('gray')
    map = flipud(colormap);
    colormap(map)
end
map = colormap;
map = [1,1,1; map];
colormap(map)
colorbar
hold off