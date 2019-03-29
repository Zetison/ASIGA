function [a,I] = sortrows2(a)

m = size(a,1);
I = uint32((1:m)');

[~,ind] = sort(a(:,2));
I = I(ind);
[~,ind] = sort(a(I,1));
I = I(ind);
clear ind
% Rearrange input rows according to the output of the sort algorithm.
a = a(I,:);


