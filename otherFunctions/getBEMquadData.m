function [Q2D_2,W2D_2,Q,W] = getBEMquadData(p_xi,p_eta,extraGP,extraGPBEM,quadMethodBEMsimpson)
p_max = max(p_xi,p_eta);
[Q2D_2,W2D_2] = tensorQuad(p_max+1+extraGPBEM,p_max+1+extraGPBEM);
W2D_2 = repmat(W2D_2,4,1); % 4 triangles around source point
if quadMethodBEMsimpson
    [Q,W] = tensorQuad(p_xi+1+extraGP,p_eta+1+extraGP);
else
    load('integration/quadData_double','quadData')
    Q = quadData.Q;
    W = quadData.W;
end