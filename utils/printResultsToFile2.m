function printResultsToFile2(varargin)
options = struct('fileEnding','.txt',...
                 'appendLabel',false,...
                 'task',[],...
                 'filename','out',...
                 'externalProvider','false',...
                 'x',[],...
                 'y',[],...
                 'xLoopName','x',...
                 'xlabel','x',...
                 'format','plain',...
                 'ylabel','y');
% if nargin > 0
%     newOptions = varargin;
%     options = updateOptions(options,newOptions);
% end
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end


x = options.x;
y = options.y;
if isempty(y)
    y = zeros(size(x,1),0);
    options.ylabel = '';
end
format = options.format;
xlabel = options.xlabel;
ylabel = options.ylabel;
xLoopName = options.xLoopName;
task = options.task;
if ~iscell(xlabel)
    xlabel = {xlabel};
end
if ~iscell(ylabel)
    ylabel = {ylabel};
end
if isempty(y)
    labels = xlabel;
else
    labels = [xlabel, ylabel];
end
if options.appendLabel
    filename = [options.filename '_' ylabel{1} 'VS' xlabel{1}];
end
filename = [options.filename, options.fileEnding];
fid = fopen(filename, 'w+','b');
switch format  
    case 'LaTeX'
        %% Create file in 'LaTeX format'
        fprintf(fid, [repmat('%%',1,30) ' Meta data ' repmat('%%',1,30) '\n']);
        fprintf(fid, '%% Boundary condition: %s\n', task.misc.BC);
        f = task.misc.omega/(2*pi);
        if numel(f) == 1
            fprintf(fid, '%% f: %0.5gHz\n', f);
        end
        fprintf(fid, '%% Incident wave direction: alpha_s = %0.5gDeg, beta_s = %0.5gDeg\n', task.ffp.alpha_s*180/pi, task.ffp.beta_s*180/pi);
        if isfield(task,'formulation')
            fprintf(fid, '%% Formulation: %s\n', task.misc.formulation);
        end
        fprintf(fid, '%% Core method: %s\n', task.misc.coreMethod);
        fprintf(fid, '%% Method: %s\n', task.misc.method);
        if ~strcmp(xLoopName,task.msh.degree)
            fprintf(fid, '%% Degree: %d\n', task.msh.degree);
        end
        if ~strcmp(xLoopName,task.msh.M)
            fprintf(fid, '%% mesh: %d\n', task.msh.M);
        end

        if isfield(task,'IEbasis')
            fprintf(fid, '%% Infinite element basis (in radial shape functions): %s\n', task.ie.IEbasis);
        end
        if isfield(task,'N') && ~strcmp(xLoopName,task.ie.N)
            fprintf(fid, '%% N: %d\n', task.ie.N);
        end
        if isfield(task.varCol{1},'totNoElems')
            fprintf(fid, '%% Total number of elements: %d\n', task.varCol{1}.totNoElems);
        end
        if isfield(task.varCol{1},'dofs')
            fprintf(fid, '%% Degrees of freedom: %d\n', task.varCol{1}.dofs);
        end
        if isfield(task.varCol{1},'timeBuildSystem')
            totTimeString = convertTimeToString(task.varCol{1}.timeBuildSystem);
            fprintf(fid, '%% Time spent building system: %s\n', totTimeString);
        end
        if isfield(task.varCol{1},'timeSolveSystem')
            totTimeString = convertTimeToString(task.varCol{1}.timeSolveSystem);
            fprintf(fid, '%% Time spent solving system: %s\n', totTimeString);
        end
        if ~isnan(xLoopName)
            fprintf(fid, '%% Basis for x-data: variation in %s\n', xLoopName);
        end

        fprintf(fid, ['%%\n' repmat('%%',1,32) ' Data ' repmat('%%',1,32) '%%\n']);
    case 'BeTSSi'
        %% Create file in 'BeTSSi format'
        if externalProvider
            fprintf(fid, '%s\n', institution); % Line 2: institution
        else
            fprintf(fid, 'NTNU_FFI\n'); % Line 2: institution
        end
        if externalProvider
            fprintf(fid, '%s\n', comments); % Line 3: Comments, Computer specs. and calculation times   
        elseif isfield(task,'actualNoDofs')
            totTimeString = convertTimeToString(task.varCol{1}.timeBuildSystem);
            fprintf(fid, 'Method: %s%s %s with %d elements (%d degrees of freedom) and NURBS degree %d. Computational time for system assembly (4x6-core Xeon 2.67 GHz): %s\n', ...
                task.misc.coreMethod, task.misc.method, task.misc.formulation, task.varCol{1}.totNoElems, task.varCol{1}.dofs, max(task.varCol{1}.nurbs{1}.degree), totTimeString); % Line 3: Comments, Computer specs. and calculation times
        else
            fprintf(fid, 'No info\n');
        end
        if isfield(task,'line4')
            fprintf(fid, '%s\n',task.line4); % Line 4: source aspect angle in degree
        else
            alpha_s = round(task.ffp.alpha_s*180/pi, 15, 'significant');
            if length(alpha_s) > 1
                fprintf(fid, '%.15g, %.15g\n',alpha_s(1),alpha_s(end)); % Line 4: source aspect angle in degree
            else
                fprintf(fid, '%.15g\n',alpha_s); % Line 4: source aspect angle in degree
            end
        end
        if isfield(task,'line5')
            fprintf(fid, '%s\n',task.line5); % Line 4: source aspect angle in degree
        else
            beta_s = round(task.ffp.beta_s*180/pi, 15, 'significant');
            if length(beta_s) > 1
                fprintf(fid, '%.15g, %.15g\n',beta_s(1),beta_s(end)); % Line 5: source elevation angle in degree
            else
                fprintf(fid, '%.15g\n',beta_s); % Line 5: source elevation angle in degree
            end
        end
        if isfield(task,'line6')
            fprintf(fid, '%s\n',task.line6); % Line 4: source aspect angle in degree
        else
            f = task.misc.omega/(2*pi);
            if strcmp(task.misc.scatteringCase, 'Sweep')
                fprintf(fid, '%.15g, %.15g\n',f(1),f(end)); % Line 6:  frequency in Hz (NOT in kHz)
            else
                fprintf(fid, '%.15g\n',f); % Line 6:  frequency in Hz (NOT in kHz)
            end
        end

        fprintf(fid, [num2str(length(x)) '\n']); % Line 7: number of lines to follow
end
fprintf(fid, [repmat('%20s',1,size(x,2)) repmat(' %20s',1,size(y,2)) '\n'], labels{:});
for i = 1:size(x,1)
    if all(isnan(x(i,:))) || (~isempty(y(i,:)) && all(isnan(y(i,:))))
        fprintf(fid,'\n');
    else
        fprintf(fid, [repmat('%20.15g',1,size(x,2)) repmat(' %20.15g',1,size(y,2)) '\n'], x(i,:), y(i,:));
    end
end
fclose(fid); 

function totTimeString = convertTimeToString(totTime)
if totTime < 60
    totTimeString = [num2str(totTime) ' seconds.'];
elseif totTime/60 < 60
    totTimeString = [num2str(totTime/60) ' minutes.'];
else
    totTimeString = [num2str(totTime/3600) ' hours.'];
end
