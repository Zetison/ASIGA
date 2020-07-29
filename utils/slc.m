function out = slc(A, ix, dim)
% out = slc(A, ix, dim) slices A in the following way
%         out = A(:, :, ..., ix, :, :,...:);
%                 ^  ^        ^  ^  
% dimensions -->  1  2       dim dim+1
if nargin < 3
    dim = 1;
end
subses = repmat({':'}, [1 ndims(A)]);
if iscell(ix)
    subses(dim) = ix;
else
    subses{dim} = ix;
end
out = A(subses{:});