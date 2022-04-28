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
% Fix invalid file names


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
        method_lc = lower(task.misc.method); % method string in lower case format
        fid = printFieldNames(fid,task,'misc');
        fid = printFieldNames(fid,task,'ffp');
        fid = printFieldNames(fid,task,'msh');
        fid = printFieldNames(fid,task,method_lc); % Only print method fields for relevant struct
        
        if isfield(task,'totNoElems') && ~strcmp(xLoopName,task.msh.M) && ~strcmp(xLoopName,task.msh.degree)
            fprintf(fid, '%% Total number of elements: %d\n', task.totNoElems);
        end
        if isfield(task,'dofs') && ~strcmp(xLoopName,task.msh.M) && ~strcmp(xLoopName,task.msh.degree)
            fprintf(fid, '%% Degrees of freedom: %d\n', task.dofs);
        end
        if isfield(task,'timeBuildSystem')
            totTimeString = convertTimeToString(task.timeBuildSystem);
            fprintf(fid, '%% Time spent building system: %s\n', totTimeString);
        end
        if isfield(task,'timeSolveSystem')
            totTimeString = convertTimeToString(task.timeSolveSystem);
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
            totTimeString = convertTimeToString(task.timeBuildSystem);
            fprintf(fid, 'Method: %s%s %s with %d elements (%d degrees of freedom) and NURBS degree %d. Computational time for system assembly (4x6-core Xeon 2.67 GHz): %s\n', ...
                task.misc.coreMethod, task.misc.method, task.misc.formulation, task.totNoElems, task.dofs, max(task.msh.degree), totTimeString); % Line 3: Comments, Computer specs. and calculation times
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

function fid = printFieldNames(fid,task,method_lc)

fieldNames = fieldnames(task.(method_lc));
for fieldname = fieldNames.'
    if ~(numel(task.(method_lc).(fieldname{1})) > 20 && ~isa(task.(method_lc).(fieldname{1}),'char'))
        switch class(task.(method_lc).(fieldname{1}))
            case 'char'
                fprintf(fid, '%% %s.%s: %s\n', method_lc, fieldname{1}, task.(method_lc).(fieldname{1}));
            case {'double','int','int32','logical'}
                fprintf(fid, '%% %s.%s: %s\n', method_lc, fieldname{1}, num2str(task.(method_lc).(fieldname{1})));
            case 'function_handle'
                fprintf(fid, '%% %s.%s: %s\n', method_lc, fieldname{1}, func2str(task.(method_lc).(fieldname{1})));
        end
    end
end
