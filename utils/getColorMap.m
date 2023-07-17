function map = getColorMap(name,ax,no_colors)

try
    RGBdata = readtable([getenv('HOME') '/kode/colormaps/' name '.rgb'],'FileType','text');
    
    map = RGBdata{:,:};
catch
    
    switch name
        case 'shsv'
            map = colormap(ax,'hsv');
            map = map(1:ceil(size(map,1)/no_colors):end,:);
            map = 0.5*(1+map); % Make minimum value 0.5
        otherwise
            error('Could not find file/not implemented')
    end
end