function f = fAnalyticPulsation_(v,C_n,y,k)

f = zeros(size(v,1),1);
for i = 1:size(y,1)
    f = f + C_n(i)/(4*pi)*exp(-1i*k*dot3(v./repmat(norm2(v),1,3),y(i,:).'));
end