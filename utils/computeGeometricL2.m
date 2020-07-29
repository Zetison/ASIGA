function [I,A] = computeGeometricL2(nurbs, newOptions)

%% Interpret input arguments

% set default values
options = struct('f', @(v)NaN(numel(v),1),...
                 'extraGP',2);

% read the acceptable names
optionNames = fieldnames(options);

if nargin == 2
    % count arguments
    nArgs = length(newOptions);
    if round(nArgs/2) ~= nArgs/2
        error('Must have propertyName/propertyValue pairs')
    end

    for pair = reshape(newOptions,2,[]) %# pair is {propName;propValue}
        inpName = pair{1}; %# make case insensitive

        if any(strcmp(inpName,optionNames))
            options.(inpName) = pair{2};
        else
            error('%s is not a recognized parameter name',inpName)
        end
    end
end

f = options.f;
extraGP = options.extraGP;

p_xi = nurbs.degree(1);
p_eta = nurbs.degree(2);

nurbspatches{1} = nurbs;
varCol.dimension = 1;
varCol = convertNURBS(nurbspatches, varCol);  
varCol = generateIGA2DMesh_new(varCol);
noElems = varCol.patches{1}.noElems;
index = varCol.patches{1}.index;
elRangeXi = varCol.patches{1}.elRange{1};
elRangeEta = varCol.patches{1}.elRange{2};
hold on
[W2D,Q2D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP);
I = 0;
A = 0;
parfor e = 1:noElems
% for e = 1:noElems
    idXi = index(e,1);
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);
    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));
    xi = parent2ParametricSpace(Xi_e, Q2D(:,1));
    eta = parent2ParametricSpace(Eta_e, Q2D(:,2));
    
    [v, dydxi, dydeta] = evaluateNURBS_2ndDeriv(nurbs, [xi eta]);

    crossProd = cross(dydxi,dydeta,2);
    J_1 = norm2(crossProd);
    fact = J_1*J_2.*W2D;
    I = I + sum(f(v).^2.*fact);
    A = A + sum(fact);
end

