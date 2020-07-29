function nurbs = read_g2(inputFileName)

fid = fopen(inputFileName,'r','b');
nurbs = cell(1,1000);
counter = 1;
while true
    A = fscanf(fid,'%d\n',4);
    if isempty(A)
        break
    end
    switch floor(A(1)/100)
        case 7
            d_p = 3;
        case 2
            d_p = 2;
        case 1
            d_p = 1;
        otherwise
            error('not implemented')
    end
        
    A = fscanf(fid,'%d\n',2);
    d = A(1);
    rational = A(2);
    number = zeros(1,d_p);
    degree = zeros(1,d_p);
    knots = cell(1,d_p);
    for i = 1:d_p
        A = fscanf(fid,'%d\n',2);
        number(i) = A(1);
        degree(i) = A(2)-1;
        knots{i} = fscanf(fid,'%f\n',number(i)+degree(i)+1)';
    end

    coeffs = zeros([d+rational,number]);
    for i = 1:prod(number)
        coeffs(:,i) = fscanf(fid,'%f\n',d+rational)';
    end
    if rational
        coeffs(1:d,:) = coeffs(1:d,:)./coeffs(end,:);
    else
        coeffs(end+1,:) = 1;
    end
    coeffs = reshape(coeffs,[d+1,number]);

    nurbs(counter) = createNURBSobject(coeffs,knots);
    counter = counter + 1;
end
fclose(fid);
nurbs(counter:end) = [];
