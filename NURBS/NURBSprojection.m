function coeffs = NURBSprojection(coeffs,d,d_p,type)

switch type
    case 'project'
        coeffs = subasgnArr(coeffs,slc(coeffs,1:d)./repmat(slc(coeffs,d+1),[d,ones(1,d_p)]),1:d);
    case 'unproject'
        coeffs = subasgnArr(coeffs,slc(coeffs,1:d).*repmat(slc(coeffs,d+1),[d,ones(1,d_p)]),1:d);
end