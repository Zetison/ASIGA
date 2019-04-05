function W = tensorWeights(Wxi,Weta,Wzeta)

switch nargin
    case 1
        W = Wxi;

    case 2
        W = Wxi*Weta.';
    case 3
        temp = Wxi*Weta.';
        W = kron(temp(:),Wzeta);
end
W = W(:);
