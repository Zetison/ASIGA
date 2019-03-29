function plotTriangulation(varCol,U,delta,xb,yb,zb,options,nameAppendix)
Lx = xb(2)-xb(1);
x = linspace(xb(1),xb(2), ceil(Lx/delta)+1);
Ly = yb(2)-yb(1);
y = linspace(yb(1),yb(2), ceil(Ly/delta)+1);
Lz = zb(2)-zb(1);
z = linspace(zb(1),zb(2), ceil(Lz/delta)+1);
min_d_Xon = 1; % minimal distance from scatterer to points in X_exterior
extraPts = 5; % extra knots in mesh for plotting on scatterer

[x,y,z] = ndgrid(x,y,z);
x = x(:);
y = y(:);
z = z(:);
DT = delaunayTriangulation(x,y,z);
X_exterior = DT.Points;

stringShift = 40;
tic
fprintf(['\n%-' num2str(stringShift) 's'], '    Triangulating the scatterer ... ')
[faces,X_on,varCol.patches,dofsToRemovePP] = triangulateNURBSsurface(varCol.patches,extraPts);
fprintf('using %12f seconds.', toc)

% tic
% fprintf(['\n%-' num2str(stringShift) 's'], '    Find points inside the scatterer ... ')
% in = intriangulation(X_on,faces,X_exterior,0);
% fprintf('using %12f seconds.', toc)

tic
fprintf(['\n%-' num2str(stringShift) 's'], '    Compute solution on the scatterer ... ')
totFieldOn = zeros(size(X_on,1), size(U,2));
noPatches = varCol.noPatches;
p_xi = varCol.degree(1); % assume p_xi is equal in all patches
p_eta = varCol.degree(2); % assume p_eta is equal in all patches

p_inc = varCol.p_inc;
omega = varCol.omega;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
noElemsPatch = varCol.noElemsPatch;

k = varCol.k;

if strcmp(varCol.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = varCol.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
end

patches = varCol.patches;
counter = 1;
for patch_i = 1:noPatches
    ppParmPts = patches{patch_i}.ppPrmPts;
    Xi = knotVecs{patch_i}{1};
    uniqueXi = unique(Xi);
    noElXi = numel(uniqueXi)-1;
    Eta = knotVecs{patch_i}{2};
    uniqueEta = unique(Eta);
    noElEta = numel(uniqueEta)-1;
    accumElems = sum(noElemsPatch(1:patch_i-1));
    
    for i = 1:size(ppParmPts,1)
        xi = ppParmPts(i,1);
        eta = ppParmPts(i,2);
        
        if xi == 1
            e_xi = noElXi;
        else
            e_xi = find(and(uniqueXi(1:end-1) <= xi, xi < uniqueXi(2:end)));
        end
        if eta == 1
            e_eta = noElEta;
        else
            e_eta = find(and(uniqueEta(1:end-1) <= eta, eta < uniqueEta(2:end)));
        end
        
        e = e_xi + (e_eta-1)*noElXi + accumElems;

        sctr = element(e,:);
        wgts = weights(element2(e,:)); % New
        pts = controlPts(sctr,:);
        U_sctr = U(sctr,:);
        
        R = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);
        y = R*pts;
        if useEnrichedBfuns
            R = R*exp(1i*k*(y*d_vec));
        end
        totFieldOn(counter,:) = R*U_sctr;
        counter = counter + 1;
    end
end
fprintf('using %12f seconds.', toc)
tic
fprintf(['\n%-' num2str(stringShift) 's'], '    Write scatterer data to file ... ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options2 = options;
options2.name = [options.name, '_surfTRI'];
options2.celltype = 'VTK_TRIANGLE';

switch varCol.method
    case 'BEM'
        scalarField = totFieldOn - p_inc(X_on);
    case 'BA'
        scalarField = totFieldOn;
end
totField = scalarField + p_inc(X_on);
data2.nodes = X_on;
data2.visElements = faces;
data2.P_inc = real(makeDynamic(p_inc(X_on), options, omega)); 
data2.scalarField = real(makeDynamic(scalarField, options, omega)); 
if options.plotAnalytic
    analytic = varCol.analytic(X_on);
    data2.analytic = real(makeDynamic(analytic, options, omega)); 
    data2.Error = abs(scalarField-analytic)./abs(analytic);
end
data2.totField = real(makeDynamic(totField, options, omega));
data2.totFieldAbs = abs(makeDynamic(totField, options, omega));

data2.omega = omega;
makeVTKfile(data2, options2);
fprintf('using %12f seconds.', toc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
fprintf(['\n%-' num2str(stringShift) 's'], '    Remove redundant points ... ')

X_on(dofsToRemovePP,:) = [];
totFieldOn(dofsToRemovePP,:) = [];

I1 = zeros(size(X_exterior,1),1,'logical');
parfor i = 1:size(X_exterior,1)
    I1(i) = any(norm2(repmat(X_exterior(i,:),size(X_on,1),1) - X_on) < delta*0.6);
end
X_exterior(I1,:) = []; 

Eps = 1e5*eps;
createMeshForCutPlanes = true;
if createMeshForCutPlanes
    I2 = all([(abs(X_exterior(:,2)) > delta+Eps), ...      % xz-plane
         	  (abs(X_exterior(:,3)) > delta+Eps), ...      % xy-plane
              (abs(X_exterior(:,1)+5) > delta+Eps), ...    % x = -5
          	  (abs(X_exterior(:,1)+18) > delta+Eps), ...   % x = -18
              (abs(X_exterior(:,1)+53) > delta+Eps)], 2);         % x = -53
    X_exterior(I2,:) = []; 
end
I3 = zeros(size(X_exterior,1),1,'logical');
parfor i = 1:size(X_exterior,1)
    I3(i) = any(norm2(repmat(X_exterior(i,:),size(X_on,1),1) - X_on) < min_d_Xon);
end
X_proximity = X_exterior(I3,:);
X_exterior = X_exterior(~I3,:);
nodes = [X_on; X_proximity; X_exterior];
DT = delaunayTriangulation(nodes);
TETR = DT.ConnectivityList;
I = zeros(size(TETR,1),1,'logical');
parfor i = 1:size(TETR,1)
    PP = nodes(TETR(i,:),:);
    I(i) = any(norm2([PP(1,:)-PP(2,:); 
                      PP(1,:)-PP(3,:); 
                      PP(1,:)-PP(4,:); 
                      PP(2,:)-PP(3,:); 
                      PP(2,:)-PP(4,:); 
                      PP(3,:)-PP(4,:)]) > 2*delta+Eps);
	
	if createMeshForCutPlanes
        normal = cross(PP(1,:)-PP(2,:),PP(1,:)-PP(3,:));
        I(i) = and(I(i),abs(dot(PP(4,:),normal) - dot(PP(3,:),normal))/dot(PP(3,:),normal) < Eps);
    end
end
TETR(I,:) = [];
fprintf('using %12f seconds.', toc)

tic
fprintf(['\n%-' num2str(stringShift) 's'], '    Compute solution in the exterior ... ')
if ~options.plotTimeOscillation
    options.plotTotFieldAbs = 1;
end

options.name = [options.name, nameAppendix];
options.celltype = 'VTK_TETRA';

data.nodes = nodes;
data.visElements = TETR;
noPtsOn = size(X_on,1);
noPtsProx = size(X_proximity,1);

switch varCol.method
    case 'BEM'
        if noPtsProx > 0
            scalarField1 = calculateScatteredPressureBEM(varCol, U, nodes(noPtsOn+1:noPtsOn+noPtsProx,:), 1, 0);
        else
            scalarField1 = [];
        end
        scalarField2 = calculateScatteredPressureBEM(varCol, U, nodes(noPtsOn+noPtsProx+1:end,:), 0, 0);
        scalarField = [totFieldOn-p_inc(X_on); scalarField1; scalarField2];
        totField = scalarField + p_inc(nodes);
    case 'BA'
        if noPtsProx > 0
            scalarField1 = calculateScatteredPressureBA(varCol, U, nodes(noPtsOn+1:noPtsOn+noPtsProx,:), 1, 0);
        else
            scalarField1 = [];
        end
        scalarField2 = calculateScatteredPressureBA(varCol, U, nodes(noPtsOn+noPtsProx+1:end,:), 0, 0);
        scalarField = [totFieldOn; scalarField1; scalarField2];
        totField = scalarField + p_inc(nodes);
%     case 'IE'
%         scalarField = calculateScatteredPressure(varCol.varColFull, varCol.U_full, nodes, 0, 0);
%         totField = scalarField + p_inc(nodes);
%     case 'IENSG'
%         scalarField = calculateScatteredPressureNonSepGeom(varCol.varColFull, varCol.U_full, nodes(noKnots+1:end,:), useExtraQuadPts, 0);
%         scalarField = [totFieldOn; scalarField];
%         totField = scalarField + p_inc(nodes);
end
fprintf('using %12f seconds.', toc)

tic
fprintf(['\n%-' num2str(stringShift) 's'], '    Write exterior solution to file ... ')
data.P_inc = real(makeDynamic(p_inc(nodes), options, omega)); 
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




