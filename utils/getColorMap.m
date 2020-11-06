function map = getColorMap(name)

RGBdata = readtable([getenv('HOME') '/kode/colormaps/' name '.rgb'],'FileType','text');

map = RGBdata{:,:};