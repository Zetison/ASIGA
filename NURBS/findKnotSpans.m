function I = findKnotSpans(degree, pt, knots)
d_p = numel(knots);
I = zeros(1,3);
for i = 1:d_p
    n = numel(knots{i})-(degree(i)+1);
    I(i) = findKnotSpan(n,degree(i), pt(i), knots{i});
end