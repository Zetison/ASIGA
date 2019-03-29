% Convert gaus quadrature weights and points from a php file to the
% GLLpoints.m file
addpath('C:\Users\Zetison\Documents\Multiprecision Computing Toolbox\')
type = 'mp';
mp.Digits(100);
n_max = 132;
[abscissa,weights] = computeGLLpoints(n_max);
for i = 2:n_max
    [temp, I] = sort(abscissa{i-1});
    temp(temp < eps) = 0;
    abscissa{i-1} = temp;
    weights{i-1} = weights{i-1}(I);
end


fid = fopen('integration\GLLpoints.m','wt+','b');
fprintf(fid,'function [points, weights] = GLLpoints(n)\n');
fprintf(fid,'points = zeros(1,n); \n');
fprintf(fid,'weights = zeros(1,n);\n');
fprintf(fid,'\n');

fprintf(fid,['if n > ' num2str(n_max) '\n']);
fprintf(fid,['\twarning(''There exist no listed gausspoints of order ' num2str(n_max+1) ' or more: using ' num2str(n_max) ' quadrature points instead...'')\n']);
fprintf(fid,['\tn = ' num2str(n_max) ';\n']);
fprintf(fid,'end\n\n');

fprintf(fid,'switch n\n');
for i = 2:n_max
    if i < 10
        fprintf(fid,'\tcase %d\n',i);
        fprintf(fid,'\t\t%% Set Gauss quadrature points\n');
        for j = floor(i/2)+1:i
            if abscissa{i-1}(j) < 0
                fprintf(fid,'\t\tpoints(%d) = %.50g;\n',j,abscissa{i-1}(j));
            else
                fprintf(fid,'\t\tpoints(%d) =  %.50g;\n',j,abscissa{i-1}(j));
            end
        end
        fprintf(fid,'\t\t\n');
        fprintf(fid,'\t\t%% Set Gauss quadrature weights\n');
        for j = floor(i/2)+1:i
            fprintf(fid,'\t\tweights(%d) = %.50g;\n',j,weights{i-1}(j));
        end
        fprintf(fid,'\t\t\n');
    else
        fprintf(fid,'\tcase %d\n',i);
        fprintf(fid,'\t\t%% Set Gauss quadrature points\n');
        for j = floor(i/2)+1:i
            if abscissa{i-1}(j) < 0
                if j < 10
                    fprintf(fid,'\t\tpoints(%d)  = %.50g;\n',j,abscissa{i-1}(j));
                else
                    fprintf(fid,'\t\tpoints(%d) = %.50g;\n',j,abscissa{i-1}(j));
                end
            else
                fprintf(fid,'\t\tpoints(%d) =  %.50g;\n',j,abscissa{i-1}(j));
            end
        end
        fprintf(fid,'\t\t\n');
        fprintf(fid,'\t\t%% Set Gauss quadrature weights\n');
        for j = floor(i/2)+1:i
            if j < 10
                fprintf(fid,'\t\tweights(%d)  = %.50g;\n',j,weights{i-1}(j));
            else
                fprintf(fid,'\t\tweights(%d) = %.50g;\n',j,weights{i-1}(j));
            end
        end
        fprintf(fid,'\t\t\n');
    end
end
fprintf(fid,'end\n\n');
fprintf(fid,'i = 1:floor(n/2);\n');
fprintf(fid,'weights(i) = weights(n-i+1);\n');
fprintf(fid,'points(i) = -points(n-i+1);');

fclose(fid);