function x = linspaceHP(a,b,n)
if n == 1
    x = a;
else
    dx = (b-a)/(n-1);
    if isa(a,'sym') || isa(b,'sym')
        x = a + vpa(0:n-1)*dx;
    elseif isa(a,'mp') || isa(b,'mp')
        x = a + mp(0:n-1)*dx;
    else
        x = a + (0:n-1)*dx;
    end
end