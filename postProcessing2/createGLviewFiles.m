function createGLviewFiles(varCol, U, extraXiPts, extraEtaPts, extraZetaPts, options, addField, field)

if nargin < 7
    addField = false;
end


Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
Zeta = varCol.nurbs.knots{3};
p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);
r = varCol.nurbs.degree(3);

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
elRangeZeta = varCol.elRangeZeta;
element = varCol.element;
weights = varCol.weights;
controlPts = varCol.controlPts;
nurbs = varCol.nurbs;
% build visualization B8 mesh

[nodes, noNodes, visElements, cornerNode, ...
          noXiKnots, noEtaKnots, noZetaKnots] = buildVisualization3dMesh_new(Xi, Eta, Zeta, extraXiPts, extraEtaPts, extraZetaPts, noElems, nurbs);

scalarField   = zeros(noNodes,1);

noUniqueXiKnots = length(unique(Xi));
noUniqueEtaKnots = length(unique(Eta));
noUniqueZetaKnots = length(unique(Zeta));

noElems = (noUniqueXiKnots-1)*(noUniqueEtaKnots-1)*(noUniqueZetaKnots-1);
for e = 1:noElems
    idXi = index(e,1);
    idEta = index(e,2);
    idZeta = index(e,3);
    
    % Modify element range such that the knotspan index do not change over
    % an element (with value 1e-12)
    
    Xi_e = elRangeXi(idXi,:); % [xi_i,xi_i+1]
    Xi_eV = linspace(Xi_e(1),Xi_e(end)-1e-12,2+extraXiPts);
    
    Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]
    Eta_eV = linspace(Eta_e(1),Eta_e(end)-1e-12,2+extraEtaPts);
    
    Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]
    Zeta_eV = linspace(Zeta_e(1),Zeta_e(end)-1e-12,2+extraZetaPts);
    
    sctr = element(e,:);          %  element scatter vector
    
    n_en = length(sctr);
    
    
    pts = controlPts(sctr,:);
        
    for idxZeta=1:2+extraZetaPts  
        for idxEta=1:2+extraEtaPts
            for idxXi = 1:2+extraXiPts            
                xi = Xi_eV(idxXi);       
                eta = Eta_eV(idxEta);  
                zeta = Zeta_eV(idxZeta);                                       
                
                [R_fun, ~, ~, ~] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
                % compute the jacobian of physical and parameter domain mapping
                % then the derivative w.r.t spatial physical coordinates
                
                nodeID = cornerNode(e) + idxXi-1   +   (idxEta-1)*noXiKnots + (idxZeta-1)*noEtaKnots*noXiKnots;
                
                if addField
                    v = R_fun*pts;
                    scalarField(nodeID) = R_fun*U(sctr) + field(v);  
                else
                    scalarField(nodeID) = R_fun*U(sctr);  
                end
 
            end
        end
    end
end

data.noElems = noElems;
data.nodes = nodes;
data.visElements = visElements;
data.scalarField = scalarField;

if exist('omega','var')
    data.omega = omega;
else
    data.omega = 1;
end

makeVTFfile_new(data, options);



