function nurbs = read_g2(inputFileName)

fid = fopen(inputFileName,'r','b');
if fid == -1
    error('Could not open file')
end
nurbs = cell(1,1000);
counter = 1;
while ~feof(fid)
    headerln = fgetl(fid);
    A = str2num(headerln);
    switch A(1)
        case 700
            d_p = 3;
        case 200
            d_p = 2;
        case 100
            d_p = 1;
        otherwise
            error('not implemented')
    end
    major = A(2);         % major version. Always 1.
    minor = A(3);         % minor version. Always 0.
    auxillary = A(4);     % auxiliary data. 0 if the default colour is chosen, 4 if the rgb-colour code is chosen.
    switch auxillary
        case 4
            color = A(5:end-1)/255;
            alpha = A(end)/255;
        case 0
            color = getColor(1);
            alpha = 1;
        otherwise
            error('Not implemented')
    end
    A = fscanf(fid,'%d\n',2);
    d = A(1);
    isRational = A(2);
    number = zeros(1,d_p);
    degree = zeros(1,d_p);
    knots = cell(1,d_p);
    for i = 1:d_p
        A = fscanf(fid,'%d\n',2);
        number(i) = A(1);
        degree(i) = A(2)-1;
        knots{i} = fscanf(fid,'%f\n',number(i)+degree(i)+1)';
    end

    coeffs = zeros([d+isRational,number]);
    for i = 1:prod(number)
        coeffs(:,i) = fscanf(fid,'%f\n',d+isRational)';
    end
    if isRational
        coeffs(1:d,:) = coeffs(1:d,:)./coeffs(end,:);
    else
        coeffs(end+1,:) = 1;
    end
    coeffs = reshape(coeffs,[d+1,number]);

    nurbs(counter) = createNURBSobject(coeffs,knots,major,minor,color,alpha);
    counter = counter + 1;
end
fclose(fid);
nurbs(counter:end) = [];
