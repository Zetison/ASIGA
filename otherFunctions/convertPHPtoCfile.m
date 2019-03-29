% Convert gaus quadrature weights and points from a php file to the
% gaussQuad.m file

weights = cell(63,1);
abscissa = cell(63,1);
fid_weights = fopen('integration\lgvalues-weights.php','r','b');
fid_abscissa = fopen('integration\lgvalues-abscissa.php','r','b');

dummy = fscanf(fid_weights,'%s\n',1)
dummy = fscanf(fid_abscissa,'%s\n',1)
dummy = fscanf(fid_weights,'%s %s %s\n',3)
dummy = fscanf(fid_abscissa,'%s %s %s\n',3)

counter = 1;
for i = 2:64
    dummy = fscanf(fid_weights,'%s %s %s\n',3);
    dummy = fscanf(fid_abscissa,'%s %s %s\n',3);
    weights{i} = zeros(1,i);
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
        weights{i-1}(j) = str2num(str_weights);
        abscissa{i-1}(j) = str2num(str_abscissa);
    end
end
fclose(fid_weights);
fclose(fid_abscissa);
for i = 2:64
    [abscissa{i-1}, I] = sort(abscissa{i-1});
    weights{i-1} = weights{i-1}(I);
end


fid = fopen('integration\gaussQuad.h','wt+','b');
fprintf(fid,'void gaussQuad(double *points, double *weights, unsigned int quadorder) {\n');
fprintf(fid,'\t// Weights with more precision may be found at\n');
fprintf(fid,'\t// http://pomax.github.io/bezierinfo/legendre-gauss.html#n2\n');
fprintf(fid,'\t\n');

fprintf(fid,'\tif(quadorder > 64){\n');
fprintf(fid,'\t\tprintf("There exist no listed Gauss points of 65 or higher");\n');
fprintf(fid,'\t\texit(EXIT_FAILURE);\n');
fprintf(fid,'\t}\n');
fprintf(fid,'\t\n');

fprintf(fid,'\tswitch(quadorder) {\n');
fprintf(fid,'\t\tcase %d:\n',1);
fprintf(fid,'\t\t\t// Set Gauss quadrature points\n');
fprintf(fid,'\t\t\tpoints[%d] = %.15g;\n',0,0);
fprintf(fid,'\t\t\t\n');
fprintf(fid,'\t\t\t// Set Gauss quadrature weights\n');
fprintf(fid,'\t\t\tweights[%d] = %.15g;\n',0,2);
fprintf(fid,'\t\t\t\n');

for i = 2:64
    if i < 10
        fprintf(fid,'\t\tcase %d:\n',i);
        fprintf(fid,'\t\t\t// Set Gauss quadrature points\n');
        for j = 1:i
            if abscissa{i-1}(j) < 0
                fprintf(fid,'\t\t\tpoints[%d] = %.15g;\n',j-1,abscissa{i-1}(j));
            else
                fprintf(fid,'\t\t\tpoints[%d] =  %.15g;\n',j-1,abscissa{i-1}(j));
            end
        end
        fprintf(fid,'\t\t\t\n');
        fprintf(fid,'\t\t\t// Set Gauss quadrature weights\n');
        for j = 1:i
            fprintf(fid,'\t\t\tweights[%d] = %.15g;\n',j-1,weights{i-1}(j));
        end
        fprintf(fid,'\t\t\t\n');
    else
        fprintf(fid,'\t\tcase %d:\n',i);
        fprintf(fid,'\t\t\t// Set Gauss quadrature points\n');
        for j = 1:i
            if abscissa{i-1}(j) < 0
                if j < 10
                    fprintf(fid,'\t\t\tpoints[%d]  = %.15g;\n',j-1,abscissa{i-1}(j));
                else
                    fprintf(fid,'\t\t\tpoints[%d] = %.15g;\n',j-1,abscissa{i-1}(j));
                end
            else
                fprintf(fid,'\t\t\tpoints[%d] =  %.15g;\n',j-1,abscissa{i-1}(j));
            end
        end
        fprintf(fid,'\t\t\t\n');
        fprintf(fid,'\t\t\t// Set Gauss quadrature weights\n');
        for j = 1:i
            if j < 10
                fprintf(fid,'\t\t\tweights[%d]  = %.15g;\n',j-1,weights{i-1}(j));
            else
                fprintf(fid,'\t\t\tweights[%d] = %.15g;\n',j-1,weights{i-1}(j));
            end
        end
        fprintf(fid,'\t\t\t\n');
    end
end
fprintf(fid,'\t}\n');
fprintf(fid,'}');

fclose(fid);