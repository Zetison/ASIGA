close all
Nxi = fluid{1}.number(1);
Neta = fluid{1}.number(2);
Nzeta = fluid{1}.number(3);
temp = reshape(U_fluid_o(1:6*Nxi*Neta*Nzeta),Nxi,Neta,Nzeta,6);
for i = 1:numel(fluid)
    N = varCol.N;
    varColPlot.plotAt = [0, 0;
                         0, 0;
                         1, 1];
    varColPlot.colorFun = @(x) log10(abs(norm2(x)));
    plotLagrange(fluid{i},[100 100 1000], 1, getColor(1), 1, temp(:,:,:,i), varColPlot);
end

axis equal
colorbar
view(60,40)