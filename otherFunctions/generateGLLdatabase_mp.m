% Convert gaus quadrature weights and points from a php file to the
% GLLpoints.m file
addpath('C:\Users\Zetison\Documents\Multiprecision Computing Toolbox\')
type = 'mp';
mp.Digits(100);
[abscissa,weights] = computeGLLpoints(64);
for i = 2:64
    [abscissa{i-1}, I] = sort(abscissa{i-1});
    weights{i-1} = weights{i-1}(I);
end


fid = fopen('integration\GLLpoints.m','wt+','b');
fprintf(fid,'function [points, weights] = GLLpoints_mp(n,useSym)\n');
fprintf(fid,'points = zeros(n,1); \n');
fprintf(fid,'weights = zeros(n,1);\n');
fprintf(fid,'if useSym\n');
fprintf(fid,'\tpoints = mp(points);\n');
fprintf(fid,'\tweights = mp(weights);\n');
fprintf(fid,'end\n');
fprintf(fid,'\n');

fprintf(fid,'if n > 64\n');
fprintf(fid,'\twarning(''There exist no listed gausspoints of order 65 or more: using 64 quadrature points instead...'')\n');
fprintf(fid,'\tn = 64;\n');
fprintf(fid,'end\n\n');

fprintf(fid,'switch n\n');
for i = 2:64
    if i < 10
        fprintf(fid,'\tcase %d\n',i);
        fprintf(fid,'\t\t%% Set Gauss quadrature points\n');
        for j = floor(i/2)+1:i
            if abscissa{i-1}(j) < 0
                fprintf(fid,'\t\tpoints(%d) = mp(''%.200g'');\n',j,abscissa{i-1}(j));
            else
                fprintf(fid,'\t\tpoints(%d) =  mp(''%.200g'');\n',j,abscissa{i-1}(j));
            end
        end
        fprintf(fid,'\t\t\n');
        fprintf(fid,'\t\t%% Set Gauss quadrature weights\n');
        for j = floor(i/2)+1:i
            fprintf(fid,'\t\tweights(%d) = mp(''%.200g'');\n',j,weights{i-1}(j));
        end
        fprintf(fid,'\t\t\n');
    else
        fprintf(fid,'\tcase %d\n',i);
        fprintf(fid,'\t\t%% Set Gauss quadrature points\n');
        for j = floor(i/2)+1:i
            if abscissa{i-1}(j) < 0
                if j < 10
                    fprintf(fid,'\t\tpoints(%d)  = mp(''%.200g'');\n',j,abscissa{i-1}(j));
                else
                    fprintf(fid,'\t\tpoints(%d) = mp(''%.200g'');\n',j,abscissa{i-1}(j));
                end
            else
                fprintf(fid,'\t\tpoints(%d) =  mp(''%.200g'');\n',j,abscissa{i-1}(j));
            end
        end
        fprintf(fid,'\t\t\n');
        fprintf(fid,'\t\t%% Set Gauss quadrature weights\n');
        for j = floor(i/2)+1:i
            if j < 10
                fprintf(fid,'\t\tweights(%d)  = mp(''%.200g'');\n',j,weights{i-1}(j));
            else
                fprintf(fid,'\t\tweights(%d) = mp(''%.200g'');\n',j,weights{i-1}(j));
            end
        end
        fprintf(fid,'\t\t\n');
    end
end
fprintf(fid,'end\n\n');
fprintf(fid,'for i = 1:floor(n/2)\n');
fprintf(fid,'\tweights(i) = weights(n-i+1);\n');
fprintf(fid,'\tpoints(i) = -points(n-i+1);\n');
fprintf(fid,'end');

fclose(fid);