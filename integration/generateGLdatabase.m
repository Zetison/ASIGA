
addpath ../../
addpath('C:\Users\Zetison\Documents\Multiprecision Computing Toolbox\')
addpath('/usr/local/AdvanpixMCT/4.6.0.13135')
startMatlabPool
mp.Digits(200);
N = 1000;
[abscissa,weights] = computeGLpoints(N);
d = mp.Digits;
Eps = 10^(-d);
for i = 1:N
    [temp, I] = sort(abscissa{i});
    temp(temp < Eps) = 0;
    x = temp;
    w = weights{i}(I);
    fid = fopen(['integration/legendreData/GL_N' num2str(i) '.csv'],'wt+','b');
    for j = floor(i/2)+1:i
        fprintf(fid,'%.100g,%.100g\n',x(j),w(j));
    end
    fclose(fid);
end
