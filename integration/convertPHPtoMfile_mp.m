% Convert gaus quadrature weights and points from a php file to the
% gaussQuad.m file
addpath('C:\Users\Zetison\Documents\Multiprecision Computing Toolbox\')
type = 'mp';
mp.Digits(400);
weights = cell(63,1);
abscissa = cell(63,1);
fid_weights = fopen('integration\lgvalues-weights.php','r','b');
fid_abscissa = fopen('integration\lgvalues-abscissa.php','r','b');

dummy = fscanf(fid_weights,'%s\n',1);
dummy = fscanf(fid_abscissa,'%s\n',1);
dummy = fscanf(fid_weights,'%s %s %s\n',3);
dummy = fscanf(fid_abscissa,'%s %s %s\n',3);

counter = 1;
for i = 2:64
    dummy = fscanf(fid_weights,'%s %s %s\n',3);
    dummy = fscanf(fid_abscissa,'%s %s %s\n',3);
    weights{i} = zeros(1,i,type);
    abscissa{i} = zeros(1,i,type);
    for j = 1:i
        str_weights = fscanf(fid_weights,'\t%s\n',1);
        str_abscissa = fscanf(fid_abscissa,'\t%s\n',1);
        if j < i
            str_weights = str_weights(1:end-1);
            str_abscissa = str_abscissa(1:end-1);
        else
            str_weights = str_weights(1:end-2);
            str_abscissa = str_abscissa(1:end-2);
        end
        weights{i-1}(j) = mp(str_weights);
        abscissa{i-1}(j) = mp(str_abscissa);
    end
end
fclose(fid_weights);
fclose(fid_abscissa);
for i = 2:64
    [abscissa{i-1}, I] = sort(abscissa{i-1});
    weights{i-1} = weights{i-1}(I);
end


fid = fopen('integration\gaussQuad_mp.m','wt+','b');
fprintf(fid,'function [points, weights] = gaussQuad_mp(n,useSym)\n');
fprintf(fid,'%% Weights with more precision may be found at\n');
fprintf(fid,'%% http://pomax.github.io/bezierinfo/legendre-gauss.html#n2\n');
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
fprintf(fid,'\tcase %d\n',1);
fprintf(fid,'\t\t%% Set Gauss quadrature points\n');
fprintf(fid,'\t\tpoints(%d) = mp(''%.200g'');\n',1,0);
fprintf(fid,'\t\t\n');
fprintf(fid,'\t\t%% Set Gauss quadrature weights\n');
fprintf(fid,'\t\tweights(%d) = mp(''%.200g'');\n',1,2);
fprintf(fid,'\t\t\n');

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