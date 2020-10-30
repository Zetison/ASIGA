function k_e = stiffnessElementMatrix(dRdX,fact,d_f,n_en,operator,C)

vectorizeOverGP = false;
switch operator
    case 'linearElasticity'
        if vectorizeOverGP % Warning this version is not generalized for anisotropic material
            idxMap = idxMapSwapCompAndBasis(d_f,n_en);
            RxRx = reshape(kron2(dRdX{1}, dRdX{1}) * fact, n_en, n_en);
            RyRy = reshape(kron2(dRdX{2}, dRdX{2}) * fact, n_en, n_en);
            RzRz = reshape(kron2(dRdX{3}, dRdX{3}) * fact, n_en, n_en);
            RyRz = reshape(kron2(dRdX{2}, dRdX{3}) * fact, n_en, n_en);
            RxRz = reshape(kron2(dRdX{1}, dRdX{3}) * fact, n_en, n_en);
            RxRy = reshape(kron2(dRdX{1}, dRdX{2}) * fact, n_en, n_en);
            mu = C(end);
            lambda = C(2,1);
            BCB_11 = (2*mu+lambda)*RxRx + mu*(RyRy+RzRz);
            BCB_22 = (2*mu+lambda)*RyRy + mu*(RxRx+RzRz);
            BCB_33 = (2*mu+lambda)*RzRz + mu*(RxRx+RyRy);
            BCB_23 = lambda*RyRz + mu*RyRz.';
            BCB_13 = lambda*RxRz + mu*RxRz.';
            BCB_12 = lambda*RxRy + mu*RxRy.';
            k_e = [BCB_11,   BCB_12,   BCB_13;
                   BCB_12.', BCB_22,   BCB_23;
                   BCB_13.', BCB_23.', BCB_33];
            k_e = k_e(idxMap); 
        else
            k_e = zeros(d_f*n_en);
            for i = 1:numel(fact)
                dRdx = [dRdX{1}(i,:); dRdX{2}(i,:); dRdX{3}(i,:)];
                B = strainDispMatrix(n_en,dRdx);
                k_e = k_e + B' * C * B * fact(i); 
            end
            temp = zeros(d_f*n_en,d_f*n_en);
            for i = 1:d_f
                for j = 1:d_f
                    temp(i:d_f:end, j:d_f:end) = k_e(1+(i-1)*n_en:i*n_en, 1+(j-1)*n_en:j*n_en);
                end
            end
            k_e = reshape(temp, (d_f*n_en)^2, 1);
        end
    case 'Laplace'
        k_e = zeros(n_en);
        if vectorizeOverGP
            for i = 1:d_f
                k_e = k_e + kron2(dRdX{i}, dRdX{i}) * fact;
            end
        else
            for i = 1:numel(fact)
                dRdx = [dRdX{1}(i,:); dRdX{2}(i,:); dRdX{3}(i,:)];
                k_e = k_e + dRdx.'*dRdx * fact(i); 
            end
            k_e = reshape(k_e, (d_f*n_en)^2, 1);
        end
    otherwise
        error('Not implemented')
end


