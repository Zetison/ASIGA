% Convert gaus quadrature weights and points from a php file to the
% gaussQuad.m file

r_arr = 0:4;
nrs = length(r_arr);
noGLpts = 50;
weights = cell(noGLpts,nrs);
abscissa = cell(noGLpts,nrs);
digits(200)
counter = 1;
for m = 1:nrs
    r = r_arr(m)
    tic
    [Q,W] = computeGaussFreudQuadratures(noGLpts,vpa(r));
    toc
    for n = 1:noGLpts
        weights{n,m} = double(W{n});
        abscissa{n,m} = double(Q{n});
    end
end
for n = 1:noGLpts
    for m = 1:nrs
        [abscissa{n,m}, I] = sort(abscissa{n,m});
        weights{n,m} = weights{n,m}(I);
    end
end


fid = fopen('integration/gaussFreudQuad.m','wt+','b');
fprintf(fid,'function [points, weights] = gaussFreudQuad(quadorder, r)\n');
fprintf(fid,'%% Weights and points calculated with the programs at\n');
fprintf(fid,'%% https://people.sc.fsu.edu/~jburkardt/m_src/gen_Freud_rule/gen_Freud_rule.html\n');
fprintf(fid,'points = zeros(quadorder,1); \n');
fprintf(fid,'weights = zeros(quadorder,1);\n');
fprintf(fid,'\n');

fprintf(fid,['if quadorder > ' num2str(noGLpts) '\n']);
fprintf(fid,['\terror(''There exist no listed Gauss-Freud points of order ' num2str(noGLpts+1) ' or more: using ' num2str(noGLpts) ' quadrature points instead...'')\n']);
fprintf(fid,['\tquadorder = ' num2str(noGLpts) ';\n']);
fprintf(fid,'end\n\n');


fprintf(fid,['if r > ' num2str(r_arr(end)) '\n']);
fprintf(fid,['\terror(''There exist no listed Gauss-Freud points of type r > ' num2str(r_arr(end)) '.'')\n']);
fprintf(fid,'end\n\n');


fprintf(fid,'switch r\n');
for m = 1:nrs
    r = r_arr(m);
    fprintf(fid,'\tcase %d\n',r);
    
    fprintf(fid,'\t\tswitch quadorder\n');

    for n = 1:noGLpts
        if n < 10
            fprintf(fid,'\t\t\tcase %d\n',n);
            fprintf(fid,'\t\t\t\t%% Set Gauss-Freud quadrature points\n');
            for j = 1:n
                fprintf(fid,'\t\t\t\tpoints(%d) = %.17g;\n',j,abscissa{n,m}(j));
            end
            fprintf(fid,'\t\t\t\t\n');
            fprintf(fid,'\t\t\t\t%% Set Gauss-Freud quadrature weights\n');
            for j = 1:n
                fprintf(fid,'\t\t\t\tweights(%d) = %.17g;\n',j,weights{n,m}(j));
            end
            fprintf(fid,'\t\t\t\t\n');
        else
            fprintf(fid,'\t\t\tcase %d\n',n);
            fprintf(fid,'\t\t\t\t%% Set Gauss-Freud quadrature points\n');
            for j = 1:n
                if j < 10
                    fprintf(fid,'\t\t\t\tpoints(%d)  = %.17g;\n',j,abscissa{n,m}(j));
                else
                    fprintf(fid,'\t\t\t\tpoints(%d) = %.17g;\n',j,abscissa{n,m}(j));
                end
            end
            fprintf(fid,'\t\t\t\t\n');
            fprintf(fid,'\t\t\t\t%% Set Gauss-Freud quadrature weights\n');
            for j = 1:n
                if j < 10
                    fprintf(fid,'\t\t\t\tweights(%d)  = %.17g;\n',j,weights{n,m}(j));
                else
                    fprintf(fid,'\t\t\t\tweights(%d) = %.17g;\n',j,weights{n,m}(j));
                end
            end
            fprintf(fid,'\t\t\t\t\n');
        end
    end
    fprintf(fid,'\t\tend\n');
end
fprintf(fid,'end');
fclose(fid);