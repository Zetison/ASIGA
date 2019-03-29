function plotOnSphere(varCol, U, options, delta, R, zb, nameAppendix)

theta_start = acos(zb(2)/R);
theta_arr = linspace(theta_start, pi-theta_start, round((pi-2*theta_start)*R/delta));
phi_arr = linspace(0,2*pi,round(2*pi*R/delta));

noXiKnots = numel(theta_arr);
noEtaKnots = numel(phi_arr);
theta = repmat(theta_arr.', noEtaKnots, 1);
phi = reshape(repmat(phi_arr, noXiKnots, 1), noXiKnots*noEtaKnots, 1);
nodes = R*[sin(theta).*cos(phi)-22.5/R, sin(theta).*sin(phi), cos(theta)];

noVisElems  = (noXiKnots-1)*(noEtaKnots-1);
visElements = zeros(noVisElems,4);
eVis = 1;

for j = 1:noEtaKnots-1
    for i = 1:noXiKnots-1
        visElements(eVis,1) = i   +   (j-1)*noXiKnots;
        visElements(eVis,2) = i+1 +   (j-1)*noXiKnots;
        visElements(eVis,3) = i+1 +       j*noXiKnots;
        visElements(eVis,4) = i   +       j*noXiKnots;

        eVis = eVis + 1;
    end
end



if ~options.plotTimeOscillation
    options.plotTotFieldAbs = 1;
end

options.name = [options.name, nameAppendix];
options.celltype = 'VTK_QUAD';

data.nodes = nodes;
data.visElements = visElements;

p_inc = varCol.p_inc;

switch varCol.method
    case 'BEM'
        scalarField = calculateScatteredPressureBEM(varCol, U, nodes, 0, 0);
        farField = calculateScatteredPressureBEM(varCol, U, nodes, 0, 1);
        totField = scalarField;
        scalarField = totField - p_inc(nodes);
    case 'BA'
        scalarField = calculateScatteredPressureBA(varCol, U, nodes, 0, 0);
        farField = calculateScatteredPressureBA(varCol, U, nodes, 0, 1);
        totField = scalarField + p_inc(nodes);
end
options.plotFarField = 1;

p_inc = varCol.p_inc;
omega = varCol.omega;

data.omega = varCol.omega;

if options.plotAnalytic
    options.plotFarFieldError = ~options.plotTimeOscillation;
    analytic = varCol.analytic(nodes);
    analyticFarField = varCol.farField(nodes);
    
    data.analytic = real(makeDynamic(analytic, options, omega)); 
    data.Error = abs(scalarField-analytic)./abs(analytic);
    data.farFieldError = abs(farField-analyticFarField)./abs(analyticFarField);
end
data.totField = real(makeDynamic(totField, options, omega));
data.totFieldAbs = abs(makeDynamic(totField, options, omega));
data.scalarField = real(makeDynamic(scalarField, options, omega)); 
data.P_inc = real(makeDynamic(p_inc(nodes), options, omega)); 
data.farField = real(makeDynamic(farField, options, omega)); 

makeVTKfile(data, options);


