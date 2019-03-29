function printResultsToFile(filename, x, y, varCol, printInBeTSSiFormat, printInLaTeXFormat, institution, comments,xlabels,ylabels)
if nargin < 9
    xlabels = {'x'};
    ylabels = {'y'};
end
if nargin > 6
    externalProvider = true;
else
    externalProvider = false;
end
if nargin < 5
    printInBeTSSiFormat = true;
end
if nargin < 6
    printInLaTeXFormat = false;
end
if printInBeTSSiFormat
    %% Create file in 'BeTSSi format'
    fid = fopen([filename '.txt'], 'w+','b');

    fprintf(fid,[varCol.saveName '.txt\n']); % Line 1: filename
    if externalProvider
        fprintf(fid, '%s\n', institution); % Line 2: institution
    else
        fprintf(fid, 'NTNU_FFI\n'); % Line 2: institution
    end
    if externalProvider
        fprintf(fid, '%s\n', comments); % Line 3: Comments, Computer specs. and calculation times   
    elseif isfield(varCol,'actualNoDofs')
        totTime = varCol.timeBuildSystem;
        if totTime < 60
            totTimeString = [num2str(totTime) ' seconds.'];
        elseif totTime/60 < 60
            totTimeString = [num2str(totTime/60) ' minutes.'];
        else
            totTimeString = [num2str(totTime/3600) ' hours.'];
        end

        fprintf(fid, 'Method: %s%s %s with %d elements (%d degrees of freedom) and NURBS degree %d. Computational time for system assembly (4x6-core Xeon 2.67 GHz): %s\n', ...
            varCol.coreMethod, varCol.method, varCol.formulation, varCol.totNoElems, varCol.actualNoDofs, max(varCol.nurbs.degree), totTimeString); % Line 3: Comments, Computer specs. and calculation times
    else
        fprintf(fid, 'No info\n');
    end

    if isfield(varCol,'line4')
        fprintf(fid, '%s\n',varCol.line4); % Line 4: source aspect angle in degree
    else
        alpha_s = round(varCol.alpha_s*180/pi, 15, 'significant');
        if length(alpha_s) > 1
            fprintf(fid, '%.15g, %.15g\n',alpha_s(1),alpha_s(end)); % Line 4: source aspect angle in degree
        else
            fprintf(fid, '%.15g\n',alpha_s); % Line 4: source aspect angle in degree
        end
    end
    if isfield(varCol,'line5')
        fprintf(fid, '%s\n',varCol.line5); % Line 4: source aspect angle in degree
    else
        beta_s = round(varCol.beta_s*180/pi, 15, 'significant');
        if length(beta_s) > 1
            fprintf(fid, '%.15g, %.15g\n',beta_s(1),beta_s(end)); % Line 5: source elevation angle in degree
        else
            fprintf(fid, '%.15g\n',beta_s); % Line 5: source elevation angle in degree
        end
    end
    if isfield(varCol,'line6')
        fprintf(fid, '%s\n',varCol.line6); % Line 4: source aspect angle in degree
    else
        if strcmp(varCol.scatteringCase, 'Sweep')
            f_arr = varCol.f_arr;
            fprintf(fid, '%.15g, %.15g\n',f_arr(1),f_arr(end)); % Line 6:  frequency in Hz (NOT in kHz)
        else
            f = varCol.f;
            fprintf(fid, '%.15g\n',f); % Line 6:  frequency in Hz (NOT in kHz)
        end
    end

    fprintf(fid, [num2str(length(x)) '\n']); % Line 7: number of lines to follow

    for i = 1:size(x,1) 
        fprintf(fid, [repmat('%20.15g,',1,size(x,2)) repmat(' %20.15g',1,size(y,2)) '\n'], x(i,:), y(i,:));
    end
    fclose(fid); 
end

if printInLaTeXFormat
    %% Create file in 'LaTeX format'
    fid = fopen([filename '.dat'], 'w+','b');
    labels = {xlabels{:}, ylabels{:}};
    fprintf(fid, [repmat('%20s',1,size(x,2)) repmat(' %20s',1,size(y,2)) '\n'], labels{:});
    
    for i = 1:size(x,1)
        fprintf(fid, [repmat('%20.15g',1,size(x,2)) repmat(' %20.15g',1,size(y,2)) '\n'], x(i,:), y(i,:));
    end
    fclose(fid); 
end