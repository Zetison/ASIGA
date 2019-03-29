function printResultsToFile2(filename, x, y, xlabel,ylabel,task,xLoopName)

if nargin < 6
    if nargin < 4
        xlabel = 'x';
        ylabel = 'y';
    end
    fid = fopen([filename '_' ylabel 'VS' xlabel '.txt'], 'w+','b');
    
    labels = [{xlabel}, {ylabel}];
    fprintf(fid, [repmat('%20s',1,size(x,2)) repmat(' %20s',1,size(y,2)) '\n'], labels{:});
    for i = 1:size(x,1)
        fprintf(fid, [repmat('%20.15g',1,size(x,2)) repmat(' %20.15g',1,size(y,2)) '\n'], x(i,:), y(i,:));
    end
    fclose(fid); 
else
    %% Create file in 'LaTeX format'
    fid = fopen([filename '_' ylabel{1} 'VS' xlabel{1} '.txt'], 'w+','b');
    labels = [xlabel, ylabel];
    fprintf(fid, [repmat('%%',1,30) ' Meta data ' repmat('%%',1,30) '\n']);
    fprintf(fid, '%% Boundary condition: %s\n', task.BC);
    if numel(task.f) == 1
        fprintf(fid, '%% f: %0.5gHz\n', task.f);
    end
    fprintf(fid, '%% Incident wave direction: alpha_s = %0.5gDeg, beta_s = %0.5gDeg\n', task.alpha_s*180/pi, task.beta_s*180/pi);
    if isfield(task,'formulation')
        fprintf(fid, '%% Formulation: %s\n', task.formulation);
    end
    fprintf(fid, '%% Core method: %s\n', task.coreMethod);
    fprintf(fid, '%% Method: %s\n', task.method);
    if nargin > 6 && ~strcmp(xLoopName,task.degree)
        fprintf(fid, '%% Degree: %d\n', task.degree);
    end
    if nargin > 6 && ~strcmp(xLoopName,task.M)
        fprintf(fid, '%% mesh: %d\n', task.M);
    end

    if isfield(task,'IEbasis')
        fprintf(fid, '%% Infinite element basis (in radial shape functions): %s\n', task.IEbasis);
    end
    if isfield(task,'N') && (nargin < 7 || ~strcmp(xLoopName,task.N))
        fprintf(fid, '%% N: %d\n', task.N);
    end
    if nargin > 6
        fprintf(fid, '%% Basis for x-data: variation in %s\n', xLoopName);
    end

    fprintf(fid, ['%%\n' repmat('%%',1,32) ' Data ' repmat('%%',1,32) '%%\n']);
    fprintf(fid, [repmat('%20s',1,size(x,2)) repmat(' %20s',1,size(y,2)) '\n'], labels{:});

    for i = 1:size(x,1)
        fprintf(fid, [repmat('%20.15g',1,size(x,2)) repmat(' %20.15g',1,size(y,2)) '\n'], x(i,:), y(i,:));
    end
    fclose(fid); 
end

