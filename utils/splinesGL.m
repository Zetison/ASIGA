function xiT = splinesGL(Xi,p)

counter = 1;

xiT = [];
uniqueXi = unique(Xi);
xiC0 = 0;
for i = 2:numel(uniqueXi)
    xi = uniqueXi(i);
    m = numel(find(Xi == xi));
    if m < p
        counter = counter + m;
    else
        GLLxi = parent2ParametricSpace([0,1],gaussLegendreQuad(p+counter).');
        if i == numel(uniqueXi)
            xiT = [xiT, xiC0+GLLxi*(xi-xiC0)];
        else
            xiT = [xiT, xiC0+GLLxi(1:end-1)*(xi-xiC0)];
        end
        counter = 1;

        xiC0 = xi;
    end
end