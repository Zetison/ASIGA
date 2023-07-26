function write_g2(nurbs,outputFileName,fmt)
if nargin < 3
    fmt = '%.16g ';
%     fmt = '%20.15g ';
    newlineChar = newline;
end
if ~(size(nurbs{1}.coeffs,1) == 1) % this is then not supplementary data
    nurbs = cleanNURBS(nurbs,[]);
end
if ~strcmp(outputFileName(end-2:end),'.g2')
    outputFileName = [outputFileName, '.g2'];
end
fid = fopen(outputFileName,'wt+','b');
for patch = 1:numel(nurbs)
    d_p = nurbs{patch}.d_p;
    controlPts = nurbs{patch}.coeffs;  
    isRational = ~all(controlPts(end,:) == 1);
    d = nurbs{patch}.d;
    switch d_p
        case 1
            classType = 100;
        case 2
            classType = 200;
        case 3
            classType = 700;
    end

    if d > 1 % assume d == 1 is the case of supplementary data
        fprintf(fid,'%d 1 0 0\n',classType);
        fprintf(fid,'%d %d\n', d, isRational);    % latter number: 1=rational, 0=non-rational

        for i = 1:d_p
            n = nurbs{patch}.number(i);
            p = nurbs{patch}.degree(i);
            Xi = nurbs{patch}.knots{i};
            fprintf(fid,'%d %d\n',n,p+1);
            temp = repmat(fmt,1,numel(Xi));
            fprintf(fid, [temp(1:end-1) '\n'], Xi);
        end

        if isRational
            for i = 1:d
                controlPts(i,:,:,:) = controlPts(i,:,:,:).*controlPts(end,:,:,:);
            end
        end
    end
    
    if isRational
        controlPts = controlPts(1:d+1,:);
    else
        controlPts = controlPts(1:d,:);
    end

    for i = 1:size(controlPts,2)
        temp = repmat(fmt,1,size(controlPts,1));
        fprintf(fid, [temp(1:end-1) '\n'], controlPts(:,i));
    end
end
fclose(fid);

