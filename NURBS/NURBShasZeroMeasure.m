function zeroMeasure = NURBShasZeroMeasure(nurbs,Eps)

if nargin < 2
    Eps = 1e-10;
end
if ~iscell(nurbs)
    nurbs = {nurbs};
end
zeroMeasure = true(size(nurbs));
for patch = 1:numel(nurbs)
    d_p = nurbs{patch}.d_p;
    d = nurbs{patch}.d;
    for iii = 1:d_p
        zeroMeasure(patch) = true;
        temp = reshape(permute(nurbs{patch}.coeffs,[1,iii+1,setdiff(1:d_p,iii)+1]),d+1,nurbs{patch}.number(iii),[]);
    
        for l = 1:size(temp,3)
            temp2 = temp(1:d,:,l).';
            uniqueCoeffs = uniquetol(temp2,Eps,'ByRows',true, 'DataScale',max(norm2(temp2)), 'OutputAllIndices', true);
            if size(uniqueCoeffs,1) ~= 1
                zeroMeasure(patch) = false;
                break
            end
        end
        if zeroMeasure(patch)
            break
        end
    end
end