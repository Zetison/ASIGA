function A = subasgnArr(A,B,idx)
if iscell(idx)
    subs = idx;
else
    subs = {[{idx},repmat({':'},1,ndims(A)-1)]};
end
A = subsasgn(A,struct('type','()','subs',subs),B);