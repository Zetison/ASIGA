function x = lineIntersection(x1,x2,x3,x4)
% Calculate intersection of two lines represented by:
% r1(t) = x1 + t*(x2-x1) and r2(s) = x3 + s*(x4-x3)
st = [x2-x1,x4-x3]\(x3-x1);
x = x1+st(1)*(x2-x1);