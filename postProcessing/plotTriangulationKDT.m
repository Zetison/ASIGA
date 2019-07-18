function plotTriangulationKDT(varCol,delta,xb,yb,zb,options,nameAppendix)
% X = triangulateRectangle([xb(1),yb(1)],[xb(2),yb(2)],ceil(1/delta)+1);
[TETR, X] = meshRectangle([xb(1),yb(1)],[xb(2),yb(2)],ceil((xb(2)-xb(1))/delta)+1);
if ~options.plotTimeOscillation
    options.plotTotFieldAbs = 1;
end

options.name = [options.name, nameAppendix];
options.celltype = 'VTK_TRIANGLE';

nodes = [X,ones(size(X,1),1)*zb(1)];
data.nodes = nodes;
data.visElements = TETR;
scalarField = calculateScatteredPressureKDT(varCol, nodes, 0);
p_INC = varCol.p_inc(nodes);
totField = scalarField + p_INC;

tic
omega = varCol.omega;
stringShift = 40;
fprintf(['\n%-' num2str(stringShift) 's'], '    Write exterior solution to file ... ')
data.P_inc = real(makeDynamic(p_INC, options, omega)); 
data.scalarField = real(makeDynamic(scalarField, options, omega)); 
if options.plotAnalytic
    analytic = varCol.analytic(nodes);
    data.analytic = real(makeDynamic(analytic, options, omega)); 
    data.Error = abs(scalarField-analytic)./abs(analytic);
end
data.totField = real(makeDynamic(totField, options, omega));
data.totFieldAbs = abs(makeDynamic(totField, options, omega));

data.omega = omega;
options.plotDisplacementVectors = false;
makeVTKfile(data, options);
fprintf('using %12f seconds.', toc)




