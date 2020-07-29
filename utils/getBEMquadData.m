function [Q2D_2,W2D_2,Q,W] = getBEMquadData(degree,extraGP,extraGPBEM,quadMethodBEM)
p_max = max(degree);
[Q2D_2, W2D_2] = gaussTensorQuad([p_max+1+extraGPBEM,p_max+1+extraGPBEM]);
W2D_2 = repmat(W2D_2,4,1); % 4 triangles around source point
switch quadMethodBEM
    case 'New'
        load('integration/quadData_double','quadData')
        Q = quadData.Q;
        W = quadData.W;
    case 'Simpson'
        [Q,W] = gaussTensorQuad(degree+1+extraGP);
    otherwise
        load('integration/quadData_double','quadData')
        Q = quadData.Q(1:(p_max+1)*2);
        W = quadData.W(1:(p_max+1)*2);
end