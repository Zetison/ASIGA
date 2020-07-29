function xiT = splinesGLL(Xi,p)

counter = 1;

xiT = [];
uniqueXi = unique(Xi);
xiC0 = 0;
for i = 2:numel(uniqueXi)
    xi = uniqueXi(i);
    if numel(find(Xi == xi)) < p
        counter = counter + 1;
    else
        GLLxi = parent2ParametricSpace([0,1],GLLpoints(p+counter));
        xiT = [xiT, xiC0+GLLxi(1:end-1)*(xi-xiC0)];
        counter = 1;
        xiC0 = xi;
    end
end
xiT = [xiT,Xi(end)];