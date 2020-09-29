function A = subasgnArr(A,B,idx,dir)
if iscell(idx)
    subs = idx;
else
    if nargin < 4
        subs = {[{idx},repmat({':'},1,ndims(A)-1)]};
    else
        subs = {repmat({':'},1,ndims(A))};
        subs{1}{dir} = idx;
    end
end
A = subsasgn(A,struct('type','()','subs',subs),B);